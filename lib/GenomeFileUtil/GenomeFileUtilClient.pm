package GenomeFileUtil::GenomeFileUtilClient;

use JSON::RPC::Client;
use POSIX;
use strict;
use Data::Dumper;
use URI;
use Bio::KBase::Exceptions;
my $get_time = sub { time, 0 };
eval {
    require Time::HiRes;
    $get_time = sub { Time::HiRes::gettimeofday() };
};

use Bio::KBase::AuthToken;

# Client version should match Impl version
# This is a Semantic Version number,
# http://semver.org
our $VERSION = "0.1.0";

=head1 NAME

GenomeFileUtil::GenomeFileUtilClient

=head1 DESCRIPTION





=cut

sub new
{
    my($class, $url, @args) = @_;
    

    my $self = {
	client => GenomeFileUtil::GenomeFileUtilClient::RpcClient->new,
	url => $url,
	headers => [],
    };

    chomp($self->{hostname} = `hostname`);
    $self->{hostname} ||= 'unknown-host';

    #
    # Set up for propagating KBRPC_TAG and KBRPC_METADATA environment variables through
    # to invoked services. If these values are not set, we create a new tag
    # and a metadata field with basic information about the invoking script.
    #
    if ($ENV{KBRPC_TAG})
    {
	$self->{kbrpc_tag} = $ENV{KBRPC_TAG};
    }
    else
    {
	my ($t, $us) = &$get_time();
	$us = sprintf("%06d", $us);
	my $ts = strftime("%Y-%m-%dT%H:%M:%S.${us}Z", gmtime $t);
	$self->{kbrpc_tag} = "C:$0:$self->{hostname}:$$:$ts";
    }
    push(@{$self->{headers}}, 'Kbrpc-Tag', $self->{kbrpc_tag});

    if ($ENV{KBRPC_METADATA})
    {
	$self->{kbrpc_metadata} = $ENV{KBRPC_METADATA};
	push(@{$self->{headers}}, 'Kbrpc-Metadata', $self->{kbrpc_metadata});
    }

    if ($ENV{KBRPC_ERROR_DEST})
    {
	$self->{kbrpc_error_dest} = $ENV{KBRPC_ERROR_DEST};
	push(@{$self->{headers}}, 'Kbrpc-Errordest', $self->{kbrpc_error_dest});
    }

    #
    # This module requires authentication.
    #
    # We create an auth token, passing through the arguments that we were (hopefully) given.

    {
	my $token = Bio::KBase::AuthToken->new(@args);
	
	if (!$token->error_message)
	{
	    $self->{token} = $token->token;
	    $self->{client}->{token} = $token->token;
	}
        else
        {
	    #
	    # All methods in this module require authentication. In this case, if we
	    # don't have a token, we can't continue.
	    #
	    die "Authentication failed: " . $token->error_message;
	}
    }

    my $ua = $self->{client}->ua;	 
    my $timeout = $ENV{CDMI_TIMEOUT} || (30 * 60);	 
    $ua->timeout($timeout);
    bless $self, $class;
    #    $self->_validate_version();
    return $self;
}




=head2 genbank_to_genome

  $result = $obj->genbank_to_genome($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a GenomeFileUtil.GenbankToGenomeParams
$result is a GenomeFileUtil.GenomeSaveResult
GenbankToGenomeParams is a reference to a hash where the following keys are defined:
	file has a value which is a GenomeFileUtil.File
	genome_name has a value which is a string
	workspace_name has a value which is a string
	source has a value which is a string
	taxon_wsname has a value which is a string
	taxon_reference has a value which is a string
	release has a value which is a string
	generate_ids_if_needed has a value which is a string
	genetic_code has a value which is an int
	type has a value which is a string
	metadata has a value which is a GenomeFileUtil.usermeta
File is a reference to a hash where the following keys are defined:
	path has a value which is a string
	shock_id has a value which is a string
	ftp_url has a value which is a string
usermeta is a reference to a hash where the key is a string and the value is a string
GenomeSaveResult is a reference to a hash where the following keys are defined:
	genome_ref has a value which is a string

</pre>

=end html

=begin text

$params is a GenomeFileUtil.GenbankToGenomeParams
$result is a GenomeFileUtil.GenomeSaveResult
GenbankToGenomeParams is a reference to a hash where the following keys are defined:
	file has a value which is a GenomeFileUtil.File
	genome_name has a value which is a string
	workspace_name has a value which is a string
	source has a value which is a string
	taxon_wsname has a value which is a string
	taxon_reference has a value which is a string
	release has a value which is a string
	generate_ids_if_needed has a value which is a string
	genetic_code has a value which is an int
	type has a value which is a string
	metadata has a value which is a GenomeFileUtil.usermeta
File is a reference to a hash where the following keys are defined:
	path has a value which is a string
	shock_id has a value which is a string
	ftp_url has a value which is a string
usermeta is a reference to a hash where the key is a string and the value is a string
GenomeSaveResult is a reference to a hash where the following keys are defined:
	genome_ref has a value which is a string


=end text

=item Description



=back

=cut

 sub genbank_to_genome
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function genbank_to_genome (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to genbank_to_genome:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'genbank_to_genome');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "GenomeFileUtil.genbank_to_genome",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'genbank_to_genome',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method genbank_to_genome",
					    status_line => $self->{client}->status_line,
					    method_name => 'genbank_to_genome',
				       );
    }
}
 


=head2 fasta_gff_to_genome

  $result = $obj->fasta_gff_to_genome($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a GenomeFileUtil.FastaGFFToGenomeParams
$result is a GenomeFileUtil.GenomeSaveResult
FastaGFFToGenomeParams is a reference to a hash where the following keys are defined:
	fasta_file has a value which is a GenomeFileUtil.File
	gff_file has a value which is a GenomeFileUtil.File
	genome_name has a value which is a string
	workspace_name has a value which is a string
	source has a value which is a string
	taxon_wsname has a value which is a string
	taxon_reference has a value which is a string
	release has a value which is a string
	genetic_code has a value which is an int
	type has a value which is a string
	metadata has a value which is a GenomeFileUtil.usermeta
File is a reference to a hash where the following keys are defined:
	path has a value which is a string
	shock_id has a value which is a string
	ftp_url has a value which is a string
usermeta is a reference to a hash where the key is a string and the value is a string
GenomeSaveResult is a reference to a hash where the following keys are defined:
	genome_ref has a value which is a string

</pre>

=end html

=begin text

$params is a GenomeFileUtil.FastaGFFToGenomeParams
$result is a GenomeFileUtil.GenomeSaveResult
FastaGFFToGenomeParams is a reference to a hash where the following keys are defined:
	fasta_file has a value which is a GenomeFileUtil.File
	gff_file has a value which is a GenomeFileUtil.File
	genome_name has a value which is a string
	workspace_name has a value which is a string
	source has a value which is a string
	taxon_wsname has a value which is a string
	taxon_reference has a value which is a string
	release has a value which is a string
	genetic_code has a value which is an int
	type has a value which is a string
	metadata has a value which is a GenomeFileUtil.usermeta
File is a reference to a hash where the following keys are defined:
	path has a value which is a string
	shock_id has a value which is a string
	ftp_url has a value which is a string
usermeta is a reference to a hash where the key is a string and the value is a string
GenomeSaveResult is a reference to a hash where the following keys are defined:
	genome_ref has a value which is a string


=end text

=item Description



=back

=cut

 sub fasta_gff_to_genome
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function fasta_gff_to_genome (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to fasta_gff_to_genome:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'fasta_gff_to_genome');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "GenomeFileUtil.fasta_gff_to_genome",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'fasta_gff_to_genome',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method fasta_gff_to_genome",
					    status_line => $self->{client}->status_line,
					    method_name => 'fasta_gff_to_genome',
				       );
    }
}
 


=head2 genome_to_gff

  $result = $obj->genome_to_gff($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a GenomeFileUtil.GenomeToGFFParams
$result is a GenomeFileUtil.GenomeToGFFResult
GenomeToGFFParams is a reference to a hash where the following keys are defined:
	genome_ref has a value which is a string
	ref_path_to_genome has a value which is a reference to a list where each element is a string
	is_gtf has a value which is a GenomeFileUtil.boolean
	target_dir has a value which is a string
boolean is an int
GenomeToGFFResult is a reference to a hash where the following keys are defined:
	gff_file has a value which is a GenomeFileUtil.File
	from_cache has a value which is a GenomeFileUtil.boolean
File is a reference to a hash where the following keys are defined:
	path has a value which is a string
	shock_id has a value which is a string
	ftp_url has a value which is a string

</pre>

=end html

=begin text

$params is a GenomeFileUtil.GenomeToGFFParams
$result is a GenomeFileUtil.GenomeToGFFResult
GenomeToGFFParams is a reference to a hash where the following keys are defined:
	genome_ref has a value which is a string
	ref_path_to_genome has a value which is a reference to a list where each element is a string
	is_gtf has a value which is a GenomeFileUtil.boolean
	target_dir has a value which is a string
boolean is an int
GenomeToGFFResult is a reference to a hash where the following keys are defined:
	gff_file has a value which is a GenomeFileUtil.File
	from_cache has a value which is a GenomeFileUtil.boolean
File is a reference to a hash where the following keys are defined:
	path has a value which is a string
	shock_id has a value which is a string
	ftp_url has a value which is a string


=end text

=item Description



=back

=cut

 sub genome_to_gff
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function genome_to_gff (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to genome_to_gff:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'genome_to_gff');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "GenomeFileUtil.genome_to_gff",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'genome_to_gff',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method genome_to_gff",
					    status_line => $self->{client}->status_line,
					    method_name => 'genome_to_gff',
				       );
    }
}
 


=head2 genome_to_genbank

  $result = $obj->genome_to_genbank($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a GenomeFileUtil.GenomeToGenbankParams
$result is a GenomeFileUtil.GenomeToGenbankResult
GenomeToGenbankParams is a reference to a hash where the following keys are defined:
	genome_ref has a value which is a string
	ref_path_to_genome has a value which is a reference to a list where each element is a string
GenomeToGenbankResult is a reference to a hash where the following keys are defined:
	genbank_file has a value which is a GenomeFileUtil.File
	from_cache has a value which is a GenomeFileUtil.boolean
File is a reference to a hash where the following keys are defined:
	path has a value which is a string
	shock_id has a value which is a string
	ftp_url has a value which is a string
boolean is an int

</pre>

=end html

=begin text

$params is a GenomeFileUtil.GenomeToGenbankParams
$result is a GenomeFileUtil.GenomeToGenbankResult
GenomeToGenbankParams is a reference to a hash where the following keys are defined:
	genome_ref has a value which is a string
	ref_path_to_genome has a value which is a reference to a list where each element is a string
GenomeToGenbankResult is a reference to a hash where the following keys are defined:
	genbank_file has a value which is a GenomeFileUtil.File
	from_cache has a value which is a GenomeFileUtil.boolean
File is a reference to a hash where the following keys are defined:
	path has a value which is a string
	shock_id has a value which is a string
	ftp_url has a value which is a string
boolean is an int


=end text

=item Description



=back

=cut

 sub genome_to_genbank
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function genome_to_genbank (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to genome_to_genbank:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'genome_to_genbank');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "GenomeFileUtil.genome_to_genbank",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'genome_to_genbank',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method genome_to_genbank",
					    status_line => $self->{client}->status_line,
					    method_name => 'genome_to_genbank',
				       );
    }
}
 


=head2 export_genome_as_genbank

  $output = $obj->export_genome_as_genbank($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a GenomeFileUtil.ExportParams
$output is a GenomeFileUtil.ExportOutput
ExportParams is a reference to a hash where the following keys are defined:
	input_ref has a value which is a string
ExportOutput is a reference to a hash where the following keys are defined:
	shock_id has a value which is a string

</pre>

=end html

=begin text

$params is a GenomeFileUtil.ExportParams
$output is a GenomeFileUtil.ExportOutput
ExportParams is a reference to a hash where the following keys are defined:
	input_ref has a value which is a string
ExportOutput is a reference to a hash where the following keys are defined:
	shock_id has a value which is a string


=end text

=item Description



=back

=cut

 sub export_genome_as_genbank
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function export_genome_as_genbank (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to export_genome_as_genbank:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'export_genome_as_genbank');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "GenomeFileUtil.export_genome_as_genbank",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'export_genome_as_genbank',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method export_genome_as_genbank",
					    status_line => $self->{client}->status_line,
					    method_name => 'export_genome_as_genbank',
				       );
    }
}
 
  
sub status
{
    my($self, @args) = @_;
    if ((my $n = @args) != 0) {
        Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
                                   "Invalid argument count for function status (received $n, expecting 0)");
    }
    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
        method => "GenomeFileUtil.status",
        params => \@args,
    });
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
                           code => $result->content->{error}->{code},
                           method_name => 'status',
                           data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
                          );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method status",
                        status_line => $self->{client}->status_line,
                        method_name => 'status',
                       );
    }
}
   

sub version {
    my ($self) = @_;
    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
        method => "GenomeFileUtil.version",
        params => [],
    });
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(
                error => $result->error_message,
                code => $result->content->{code},
                method_name => 'export_genome_as_genbank',
            );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(
            error => "Error invoking method export_genome_as_genbank",
            status_line => $self->{client}->status_line,
            method_name => 'export_genome_as_genbank',
        );
    }
}

sub _validate_version {
    my ($self) = @_;
    my $svr_version = $self->version();
    my $client_version = $VERSION;
    my ($cMajor, $cMinor) = split(/\./, $client_version);
    my ($sMajor, $sMinor) = split(/\./, $svr_version);
    if ($sMajor != $cMajor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Major version numbers differ.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor < $cMinor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Client minor version greater than Server minor version.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor > $cMinor) {
        warn "New client version available for GenomeFileUtil::GenomeFileUtilClient\n";
    }
    if ($sMajor == 0) {
        warn "GenomeFileUtil::GenomeFileUtilClient version is $svr_version. API subject to change.\n";
    }
}

=head1 TYPES



=head2 boolean

=over 4



=item Description

A boolean - 0 for false, 1 for true.
@range (0, 1)


=item Definition

=begin html

<pre>
an int
</pre>

=end html

=begin text

an int

=end text

=back



=head2 File

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
path has a value which is a string
shock_id has a value which is a string
ftp_url has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
path has a value which is a string
shock_id has a value which is a string
ftp_url has a value which is a string


=end text

=back



=head2 usermeta

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the key is a string and the value is a string
</pre>

=end html

=begin text

a reference to a hash where the key is a string and the value is a string

=end text

=back



=head2 GenbankToGenomeParams

=over 4



=item Description

genome_name - becomes the name of the object
workspace_name - the name of the workspace it gets saved to.
source - Source of the file typically something like RefSeq or Ensembl
taxon_ws_name - where the reference taxons are : ReferenceTaxons
    taxon_reference - if defined, will try to link the Genome to the specified
taxonomy object insteas of performing the lookup during upload
release - Release or version number of the data 
  per example Ensembl has numbered releases of all their data: Release 31
generate_ids_if_needed - If field used for feature id is not there, 
  generate ids (default behavior is raising an exception)
genetic_code - Genetic code of organism. Overwrites determined GC from 
  taxon object
type - Reference, Representative or User upload


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
file has a value which is a GenomeFileUtil.File
genome_name has a value which is a string
workspace_name has a value which is a string
source has a value which is a string
taxon_wsname has a value which is a string
taxon_reference has a value which is a string
release has a value which is a string
generate_ids_if_needed has a value which is a string
genetic_code has a value which is an int
type has a value which is a string
metadata has a value which is a GenomeFileUtil.usermeta

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
file has a value which is a GenomeFileUtil.File
genome_name has a value which is a string
workspace_name has a value which is a string
source has a value which is a string
taxon_wsname has a value which is a string
taxon_reference has a value which is a string
release has a value which is a string
generate_ids_if_needed has a value which is a string
genetic_code has a value which is an int
type has a value which is a string
metadata has a value which is a GenomeFileUtil.usermeta


=end text

=back



=head2 GenomeSaveResult

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
genome_ref has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
genome_ref has a value which is a string


=end text

=back



=head2 FastaGFFToGenomeParams

=over 4



=item Description

genome_name - becomes the name of the object
workspace_name - the name of the workspace it gets saved to.
source - Source of the file typically something like RefSeq or Ensembl
taxon_ws_name - where the reference taxons are : ReferenceTaxons
taxon_reference - if defined, will try to link the Genome to the specified
taxonomy object insteas of performing the lookup during upload
release - Release or version number of the data 
  per example Ensembl has numbered releases of all their data: Release 31
generate_ids_if_needed - If field used for feature id is not there, 
  generate ids (default behavior is raising an exception)
genetic_code - Genetic code of organism. Overwrites determined GC from 
  taxon object
type - Reference, Representative or User upload


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
fasta_file has a value which is a GenomeFileUtil.File
gff_file has a value which is a GenomeFileUtil.File
genome_name has a value which is a string
workspace_name has a value which is a string
source has a value which is a string
taxon_wsname has a value which is a string
taxon_reference has a value which is a string
release has a value which is a string
genetic_code has a value which is an int
type has a value which is a string
metadata has a value which is a GenomeFileUtil.usermeta

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
fasta_file has a value which is a GenomeFileUtil.File
gff_file has a value which is a GenomeFileUtil.File
genome_name has a value which is a string
workspace_name has a value which is a string
source has a value which is a string
taxon_wsname has a value which is a string
taxon_reference has a value which is a string
release has a value which is a string
genetic_code has a value which is an int
type has a value which is a string
metadata has a value which is a GenomeFileUtil.usermeta


=end text

=back



=head2 GenomeToGFFParams

=over 4



=item Description

is_gtf - optional flag switching export to GTF format (default is 0, 
    which means GFF)
target_dir - optional target directory to create file in (default is
    temporary folder with name 'gff_<timestamp>' created in scratch)


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
genome_ref has a value which is a string
ref_path_to_genome has a value which is a reference to a list where each element is a string
is_gtf has a value which is a GenomeFileUtil.boolean
target_dir has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
genome_ref has a value which is a string
ref_path_to_genome has a value which is a reference to a list where each element is a string
is_gtf has a value which is a GenomeFileUtil.boolean
target_dir has a value which is a string


=end text

=back



=head2 GenomeToGFFResult

=over 4



=item Description

from_cache is 1 if the file already exists and was just returned, 0 if
the file was generated during this call.


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
gff_file has a value which is a GenomeFileUtil.File
from_cache has a value which is a GenomeFileUtil.boolean

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
gff_file has a value which is a GenomeFileUtil.File
from_cache has a value which is a GenomeFileUtil.boolean


=end text

=back



=head2 GenomeToGenbankParams

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
genome_ref has a value which is a string
ref_path_to_genome has a value which is a reference to a list where each element is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
genome_ref has a value which is a string
ref_path_to_genome has a value which is a reference to a list where each element is a string


=end text

=back



=head2 GenomeToGenbankResult

=over 4



=item Description

from_cache is 1 if the file already exists and was just returned, 0 if
the file was generated during this call.


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
genbank_file has a value which is a GenomeFileUtil.File
from_cache has a value which is a GenomeFileUtil.boolean

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
genbank_file has a value which is a GenomeFileUtil.File
from_cache has a value which is a GenomeFileUtil.boolean


=end text

=back



=head2 ExportParams

=over 4



=item Description

input and output structure functions for standard downloaders


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
input_ref has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
input_ref has a value which is a string


=end text

=back



=head2 ExportOutput

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
shock_id has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
shock_id has a value which is a string


=end text

=back



=cut

package GenomeFileUtil::GenomeFileUtilClient::RpcClient;
use base 'JSON::RPC::Client';
use POSIX;
use strict;

#
# Override JSON::RPC::Client::call because it doesn't handle error returns properly.
#

sub call {
    my ($self, $uri, $headers, $obj) = @_;
    my $result;


    {
	if ($uri =~ /\?/) {
	    $result = $self->_get($uri);
	}
	else {
	    Carp::croak "not hashref." unless (ref $obj eq 'HASH');
	    $result = $self->_post($uri, $headers, $obj);
	}

    }

    my $service = $obj->{method} =~ /^system\./ if ( $obj );

    $self->status_line($result->status_line);

    if ($result->is_success) {

        return unless($result->content); # notification?

        if ($service) {
            return JSON::RPC::ServiceObject->new($result, $self->json);
        }

        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    elsif ($result->content_type eq 'application/json')
    {
        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    else {
        return;
    }
}


sub _post {
    my ($self, $uri, $headers, $obj) = @_;
    my $json = $self->json;

    $obj->{version} ||= $self->{version} || '1.1';

    if ($obj->{version} eq '1.0') {
        delete $obj->{version};
        if (exists $obj->{id}) {
            $self->id($obj->{id}) if ($obj->{id}); # if undef, it is notification.
        }
        else {
            $obj->{id} = $self->id || ($self->id('JSON::RPC::Client'));
        }
    }
    else {
        # $obj->{id} = $self->id if (defined $self->id);
	# Assign a random number to the id if one hasn't been set
	$obj->{id} = (defined $self->id) ? $self->id : substr(rand(),2);
    }

    my $content = $json->encode($obj);

    $self->ua->post(
        $uri,
        Content_Type   => $self->{content_type},
        Content        => $content,
        Accept         => 'application/json',
	@$headers,
	($self->{token} ? (Authorization => $self->{token}) : ()),
    );
}



1;
