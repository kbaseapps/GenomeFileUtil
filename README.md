# GenomeFileUtil
This is the basic readme for this module. Include any usage or deployment instructions and links to other documentation here.

## Current Status

| Branch  | Build                                                              | Coverage                                                                         | LGTM Alerts                                                     |
| ------- | ------------------------------------------------------------------ | -------------------------------------------------------------------------------- | --------------------------------------------------------------- |
| master  | [![KBase SDK Tests](https://github.com/kbaseapps/GenomeFileUtil/workflows/KBase%20SDK%20Tests/badge.svg)](https://github.com/kbaseapps/GenomeFileUtil/actions?query=workflow%3A%22KBase+SDK+Tests%22)  | [![codecov](https://codecov.io/gh/kbaseapps/GenomeFileUtil/branch/master/graph/badge.svg)](https://codecov.io/gh/kbaseapps/GenomeFileUtil)  | [![Total alerts](https://img.shields.io/lgtm/alerts/g/kbaseapps/GenomeFileUtil.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/kbaseapps/GenomeFileUtil/alerts/)  |


## Details

### Testing

The tests use data from KBase's CI environment. In `test_local/test.cfg`, set the following:

```
test_token={ci_dev_token}
# a second test token from a different account is required
test_token2={ci_dev_token2}
kbase_endpoint=https://ci.kbase.us/services
auth_service_url=https://ci.kbase.us/services/auth/api/legacy/KBase/Sessions/Login
```
