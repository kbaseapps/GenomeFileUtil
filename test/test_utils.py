from copy import deepcopy
import hashlib
import json
import os
import re
from datetime import datetime
import uuid

_UPA_PATTERN = re.compile(r'^[0-9]+\/[0-9]+\/[0-9]$')
_OBJECT_VERSION_PATTERN = re.compile(r'^[0-9]+\/[0-9]+\/1$')
_PROV_SUBACTION_VERSION_PATTERN = re.compile(r'^\d+\.\d+\.\d+-(release|beta)$')

PROVENANCE = [
    {
        "service": "GenomeFileUtil",
        "service_ver": "local-dev",
        "method": "run_local_tests",
        "method_params": [],
        "input_ws_objects": [],
        "resolved_ws_objects": [],
        "intermediate_incoming": [],
        "intermediate_outgoing": [],
        "external_data": [],
        "subactions": [
            {
                "name": "GenomeFileUtil",
                "ver": "local-dev",
                "code_url": "https://localhost",
            },
            {
                "name": "AssemblyUtil",
                "code_url": "https://github.com/kbaseapps/AssemblyUtil",
            },
            {
                "name": "DataFileUtil",
                "code_url": "https://github.com/kbaseapps/DataFileUtil",
            },
            # {
            #     "name": "GenomeAnnotationAPI",
            #     "code_url": "https://github.com/kbase/genome_annotation_api",
            # },
            {
                "name": "WsLargeDataIO",
                "code_url": "https://github.com/kbaseapps/WsLargeDataIO",
            },
        ],
        "custom": {},
        "description": "KBase SDK method run via the KBase Execution Engine",
    }
]

def check_result_object_info_provenance_data(
    results,
    file_names,
    scratch_dir,
    ws_client,
    hs_client,
    dfu_client,
    expected_wsID,
    expected_wsName,
    expected_metadata,
    expected_provenance,
    expected_data,
    expected_md5sum,
    is_genome=True,
):
    for idx, res in enumerate(results):
        expected_info, ref, expected_assembly_upa = _get_info_and_ref(res, is_genome)
        obj = _get_object(ws_client, ref)
        _check_info(obj, file_names[idx], expected_metadata[idx], expected_wsID, expected_wsName, expected_info, expected_assembly_upa, is_genome)
        _check_prov(obj, expected_provenance)
        _check_data(obj, file_names[idx], scratch_dir, hs_client, dfu_client, expected_data[idx], expected_md5sum[idx], expected_assembly_upa, is_genome)

def _get_info_and_ref(result, is_genome):
    # Process the result returned by import_genbank_mass
    if "genome_ref" in result and "assembly_ref" in result:
        ref_type = "genome_ref" if is_genome else "assembly_ref"
        object_info = "genome_info" if is_genome else "assembly_info"

        info = result[object_info]
        ref = result[ref_type]
        assembly_upa = result["assembly_ref"]

    # Process the result returned by save_genome_mass
    elif "info" in result:
        info = result["info"]
        ref = f'{info[6]}/{info[0]}/{info[4]}'
        assembly_upa = info[10]["Assembly Object"]

    else:
        raise ValueError("Result must originate from either import_genbank_mass or save_genome_mass.")

    return info, ref, assembly_upa

def _get_object(ws_client, ref):
    assert _OBJECT_VERSION_PATTERN.match(ref)
    return ws_client.get_objects2({"objects": [{'ref': ref}]})["data"][0]

def _check_info(obj, file_name, expected_metadata, expected_wsID, expected_wsName, expected_info, expected_assembly_upa, is_genome):
    info = obj["info"]
    object_name = file_name if is_genome else file_name + "_assembly"
    object_type = 'KBaseGenomes.Genome' if is_genome else 'KBaseGenomeAnnotations.Assembly'
    retrieved_metadata = _retrieve_genome_metadata(info[10], expected_assembly_upa) if is_genome else info[10]

    assert info == expected_info
    assert info[1] == object_name
    assert info[2].split('-')[0] == object_type

    # check version
    assert info[4] == 1
    # datetime.fromisoformat is not available in Python 3.6 or below
    assert datetime.strptime(info[3], '%Y-%m-%dT%H:%M:%S+%f')
    assert info[6] == expected_wsID
    assert info[7] == expected_wsName

    # check metadata
    assert retrieved_metadata == expected_metadata

def _retrieve_genome_metadata(metadata, expected_assembly_upa):
    # make a deep copy to avoid modifying the original metadata
    metadata = deepcopy(metadata)
    retrieved_assembly_upa = metadata.pop("Assembly Object")
    _check_assembly_upa(retrieved_assembly_upa, expected_assembly_upa)
    return metadata

def _check_assembly_upa(retrieved_upa, expected_upa):
    assert _UPA_PATTERN.match(retrieved_upa)
    assert retrieved_upa == expected_upa

def _check_prov(obj, expected_provenance):
    provenance = obj["provenance"]
    retrieved_provenance = _retrieve_provenance(provenance)

    print(f"expected_provenance is: {expected_provenance}")
    print("*" * 30)
    print(f"retrieved_provenance is: {retrieved_provenance}")

    assert retrieved_provenance == expected_provenance

def _retrieve_provenance(provenance):
    # make a deep copy to avoid modifying the original provenance
    provs = [prov.copy() for prov in provenance]
    for prov in provs:
        for key in ["time", "epoch"]:
            prov.pop(key)
        for subaction in prov['subactions']:
            subaction.pop("commit")
            if subaction["name"] != "GenomeFileUtil":
                version = subaction.pop("ver")
                assert _PROV_SUBACTION_VERSION_PATTERN.match(version)
    return provs

def _check_data(
    obj,
    file_name,
    scratch_dir,
    hs_client,
    dfu_client,
    expected_data,
    expected_md5sum,
    expected_assembly_ref,
    is_genome
):
    data = obj["data"]

    retrieved_data, retrieved_md5sum, retrieved_node_filename = (
        _retrieve_genome_data(dfu_client, scratch_dir, data, expected_assembly_ref)
        if is_genome
        else _retrieve_assembly_data(hs_client, dfu_client, scratch_dir, data)
    )

    expected_node_filename = (
        file_name + ".gz"
        if is_genome
        else file_name + "_assembly.fasta"
    )

    assert retrieved_data == expected_data
    assert retrieved_md5sum == expected_md5sum
    assert retrieved_node_filename == expected_node_filename

def _retrieve_genome_data(dfu_client, scratch_dir, data, expected_assembly_ref):
    # make a deep copy to avoid modifying the original genome data
    data = deepcopy(data)

    for key in ["cdss", "features", "mrnas", "non_coding_features"]:
        for dist in data.get(key):
            if dist.get("aliases"):
                dist["aliases"] = sorted(dist["aliases"])

    retrieved_assembly_ref = data.pop("assembly_ref")
    assert _UPA_PATTERN.match(retrieved_assembly_ref)
    assert retrieved_assembly_ref == expected_assembly_ref

    # check handle ref
    handle_id = data.pop("genbank_handle_ref")
    file_ret = _download_file_from_blobstore(dfu_client, scratch_dir, handle_id)

    file_path = file_ret['file_path']
    retrieved_node_filename = file_ret['node_file_name']
    retrieved_genome_md5sum = calculate_md5sum(file_path)

    for ontology_event in data.get("ontology_events", []):
        ontology_event.pop("timestamp")
        ontology_ref = ontology_event.pop("ontology_ref")
        assert _UPA_PATTERN.match(ontology_ref)

    return data, retrieved_genome_md5sum, retrieved_node_filename

def _retrieve_assembly_data(hs_client, dfu_client, scratch_dir, data):
    # make a deep copy to avoid modifying the original assembly data
    data = deepcopy(data)

    handle_id = data.pop("fasta_handle_ref")
    assert handle_id.split("_")[0] == "KBH"

    handle_info = data["fasta_handle_info"]
    blob_id = handle_info.pop("shock_id")

    handle = handle_info['handle']
    assert handle.pop('hid') == handle_id
    assert handle.pop('id') == blob_id

    # check handle_id and blob_id are match
    retrieved_blob_id = _get_blob_id(hs_client, handle_id)
    assert retrieved_blob_id == blob_id

    # check handle ref
    file_ret = _download_file_from_blobstore(dfu_client, scratch_dir, handle_id)
    file_path = file_ret['file_path']
    retrieved_node_filename = file_ret['node_file_name']
    retrieved_assembly_md5sum = calculate_md5sum(file_path)
    assert retrieved_assembly_md5sum == handle["remote_md5"]

    url = handle.get('url')
    assert url.startswith('https://')
    assert url.endswith('kbase.us/services/shock-api')

    return data, retrieved_assembly_md5sum, retrieved_node_filename

def _download_file_from_blobstore(dfu_client, scratch_dir, handle_id):
    output_dir = scratch_dir + "/" + str(uuid.uuid4())
    os.makedirs(output_dir)
    file_ret = dfu_client.shock_to_file(
        {
            'handle_id': handle_id,
            'file_path': output_dir,
            'unpack': 'unpack',
        }
    )
    return file_ret

def calculate_md5sum(file_path):
    md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            md5.update(chunk)
    return md5.hexdigest()

def load_expected_data(json_path):
    with open(json_path, "r") as read_file:
        data = json.load(read_file)
    return data

def _get_blob_id(hs_client, handle_id):
    handles = hs_client.hids_to_handles([handle_id])
    blob_id = handles[0]['id']
    return blob_id
