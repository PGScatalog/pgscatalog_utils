import pytest
import numpy as np

from pgscatalog_utils.validate.formatted.validator import init_validator as formatted_init_validator
from pgscatalog_utils.validate.harmonized_position.validator import init_validator as hmpos_init_validator


log_file = 'VALIDATE.log'


###### Formatted scoring files ######
def _get_formatted_validator(test_file):
    validator = formatted_init_validator(test_file,log_file,None)
    return validator

def _valid_file(test_file):
    validator = _get_formatted_validator(test_file)
    assert validator.validate_file_extension()
    assert validator.validate_headers()
    assert validator.validate_data()
    assert validator.is_file_valid()

def _failed_file(test_file):
    validator = _get_formatted_validator(test_file)
    assert validator.validate_file_extension()
    assert validator.validate_headers()
    assert not validator.validate_data()
    assert not validator.is_file_valid()

def _failed_header_file(test_file):
    validator = _get_formatted_validator(test_file)
    assert validator.validate_file_extension()
    validator.header = np.delete(validator.header,np.s_[0,1,2])
    assert not validator.validate_headers()


# Valid file with rsID, chr_name and chr_position
def test_valid_formatted_file_rsID_and_pos(test_file_1):
    _valid_file(test_file_1)

# Valid file with rsID only
def test_valid_formatted_file_rsID_only(test_file_2):
    _valid_file(test_file_2)

# Valid file with chr_name and chr_position
def test_valid_formatted_file_pos_only(test_file_3):
    _valid_file(test_file_3)

# File made invalid file by removing some mandatory column headers
def test_failed_formatted_file_missing_header(test_file_1):
    _failed_header_file(test_file_1)

# Invalid file with several data content issues
def test_failed_formatted_file_data_issues(test_file_4):
    _failed_file(test_file_4)



###### Harmonized (Position) scoring files ######
def _get_hmpos_validator(test_file):
    validator = hmpos_init_validator(test_file,log_file,None)
    return validator

def _valid_hmpos_file(test_file):
    validator = _get_hmpos_validator(test_file)
    assert validator.validate_file_extension()
    assert validator.validate_headers()
    assert validator.validate_data()
    assert validator.is_file_valid()

def _failed_file(test_file):
    validator = _get_formatted_validator(test_file)
    assert validator.validate_file_extension()
    assert validator.validate_headers()
    assert not validator.validate_data()
    assert not validator.is_file_valid()

def _failed_header_file(test_file):
    validator = _get_formatted_validator(test_file)
    assert validator.validate_file_extension()
    validator.header = np.delete(validator.header,np.s_[0,1,2])
    assert not validator.validate_headers()


## GRCh37 ## 
# Valid file with rsID, chr_name and chr_position
def test_valid_hmpos_file_rsID_and_pos_37(test_hmpos_file_GRCh37_1):
    _valid_hmpos_file(test_hmpos_file_GRCh37_1)
# Valid file with rsID only
def test_valid_formatted_file_rsID_only_37(test_hmpos_file_GRCh37_2):
    _valid_hmpos_file(test_hmpos_file_GRCh37_2)
# Valid file with chr_name and chr_position
def test_valid_formatted_file_pos_only_37(test_hmpos_file_GRCh37_3):
    _valid_file(test_hmpos_file_GRCh37_3)

## GRCh38 ##
# Valid file with rsID, chr_name and chr_position
def test_valid_hmpos_file_rsID_and_pos_38(test_hmpos_file_GRCh38_1):
    _valid_hmpos_file(test_hmpos_file_GRCh38_1)
# Valid file with rsID only
def test_valid_formatted_file_rsID_only_38(test_hmpos_file_GRCh38_2):
    _valid_hmpos_file(test_hmpos_file_GRCh38_2)
# Valid file with chr_name and chr_position
def test_valid_formatted_file_pos_only_38(test_hmpos_file_GRCh38_3):
    _valid_file(test_hmpos_file_GRCh38_3)


######################################################

@pytest.fixture
def test_file_1():
   return './data/test_scoring_file_1.txt.gz'

@pytest.fixture
def test_file_2():
   return './data/test_scoring_file_2.txt.gz'

@pytest.fixture
def test_file_3():
   return './data/test_scoring_file_3.txt.gz'

@pytest.fixture
def test_file_4():
   return './data/test_scoring_file_4.txt.gz'

@pytest.fixture
def test_hmpos_file_GRCh37_1():
   return './data/test_scoring_file_hmpos_37_1.txt.gz'

@pytest.fixture
def test_hmpos_file_GRCh38_1():
   return './data/test_scoring_file_hmpos_38_1.txt.gz'

@pytest.fixture
def test_hmpos_file_GRCh37_2():
   return './data/test_scoring_file_hmpos_37_2.txt.gz'

@pytest.fixture
def test_hmpos_file_GRCh38_2():
   return './data/test_scoring_file_hmpos_38_2.txt.gz'

@pytest.fixture
def test_hmpos_file_GRCh37_3():
   return './data/test_scoring_file_hmpos_37_3.txt.gz'

@pytest.fixture
def test_hmpos_file_GRCh38_3():
   return './data/test_scoring_file_hmpos_38_3.txt.gz'