import re
from pgscatalog_utils.validate.schemas import *
from pgscatalog_utils.validate.validator_base import *

'''
PGS Catalog Harmonized file validator
- using pandas_schema https://github.com/TMiguelT/PandasSchema
'''

class ValidatorPos(ValidatorBase):
    ''' Validator for the HmPOS Harmonized file format. '''

    def __init__(self, file, score_dir=None, logfile="VALIDATE.log", error_limit=0):
        super().__init__(file, score_dir, logfile, error_limit)
        self.meta_format = HM_META_POS
        self.validators = POS_VALIDATORS
        self.valid_cols = VALID_COLS_POS
        self.valid_type = VALID_TYPE_POS


    def extract_specific_metadata(self,line):
        ''' Extract some of the metadata. '''
        match_variants_number = re.search(r'#variants_number=(\d+)', line)
        if match_variants_number:
            self.variants_number = int(match_variants_number.group(1))


    def validate_line_content(self,cols_content,var_line_number):
        ''' Populate the abstract method from ValidatorBase, to check some data in esch row. '''
        # Check lines
        line_dict = dict(zip(self.header, cols_content))
        line_cols = line_dict.keys()
        # Check each chromosome data is consistent
        chr_cols = ['chr_name', 'hm_chr', 'hm_match_chr']
        if all(col_name in line_cols for col_name in chr_cols):
            if line_dict['chr_name'] == line_dict['hm_chr'] and line_dict['hm_match_chr'] != 'True':
                self.logger.error(f"- Variant line {var_line_number} | 'hm_match_chr' should be 'True': same chromosome ('chr_name={line_dict['chr_name']}' vs 'hm_chr={line_dict['hm_chr']}')")
        # Check each position data is consistent
        pos_cols = ['chr_position', 'hm_pos', 'hm_match_pos']
        if all(col_name in line_cols for col_name in pos_cols):
            if line_dict['chr_position'] == line_dict['hm_pos'] and line_dict['hm_match_pos'] != 'True':
                self.logger.error(f"- Variant line {var_line_number} | 'hm_match_pos' should be 'True': same position ('chr_position={line_dict['chr_position']}' vs 'hm_pos={line_dict['hm_pos']}')")


    def validate_filename(self):
        ''' Validate the file name structure. '''
        pgs_id, build = None, None
        # hmPOS
        filename = self.file.split('/')[-1].split('.')[0]
        filename_parts = filename.split('_hmPOS_')
        if len(filename_parts) != 2:
            self.logger.error("Filename: {} should follow the pattern <pgs_id>_hmPOS_<build>.txt.gz [build=GRChXX]".format(filename))
            return False
        else:
            pgs_id, build = filename_parts
        self.file_pgs_id = pgs_id
        self.file_genomebuild = build
        if not self.check_build_is_legit(build):
            self.logger.error("Build: {} is not an accepted build value".format(build))
            return False
        self.logger.info("Filename looks good!")
        return True


    def validate_headers(self):
        ''' Validate the list of column names. '''
        # Check if it has at least a "SNP" column or a "chromosome" column
        self.setup_field_validation()
        required_is_subset = set(STD_COLS_VAR_POS).issubset(self.header)
        if not required_is_subset:
            self.logger.error("Required headers: {} are not in the file header: {}".format(STD_COLS_VAR_POS, self.header))
       
        # Check if it has at least a "SNP" column or a "chromosome" column
        required_pos = set(SNP_COLS_VAR_POS).issubset(self.header)
        if not required_pos:
            # check if everything but snp:
            required_pos = set(CHR_COLS_VAR_POS).issubset(self.header)
            if not required_pos:
                self.logger.error("One of the following required header is missing: '{}' and/or '{}' are not in the file header: {}".format(SNP_COLS_VAR_POS, CHR_COLS_VAR_POS, self.header))
                required_is_subset = required_pos
        
        return required_is_subset


##################################################################

def init_validator(file, logfile, score_dir=None) -> ValidatorPos:
    validator = ValidatorPos(file=file, score_dir=score_dir, logfile=logfile)
    return validator

# def run_validator(file, check_filename, logfile, score_dir=None):

#     validator = ValidatorPos(file=file, score_dir=score_dir, logfile=logfile)

#     validator.logger.propagate = False

#     if not file or not logfile:
#         validator.logger.info("Missing file and/or logfile")
#         validator.logger.info("Exiting before any further checks")
#         sys.exit()
#     if not os.path.exists(file):
#         validator.logger.info("Error: the file '"+file+"' can't be found")
#         validator.logger.info("Exiting before any further checks")
#         sys.exit()

#     is_ok_to_run_validation = 1
#     validator.logger.info("Validating file extension...")
#     if not validator.validate_file_extension():
#         validator.logger.info("Invalid file extension: {}".format(file))
#         validator.logger.info("Exiting before any further checks")
#         is_ok_to_run_validation = 0

#     if is_ok_to_run_validation and check_filename:
#         validator.logger.info("Validating file name...")
#         if not validator.validate_filename():
#             validator.logger.info("Invalid filename: {}".format(file))
#             is_ok_to_run_validation = 0

#     if is_ok_to_run_validation:
#         validator.logger.info("Comparing filename with metadata...")
#         if not validator.compare_with_filename():
#             validator.logger.info("Discrepancies between filename information and metadata: {}".format(file))
#             is_ok_to_run_validation = 0

#     if is_ok_to_run_validation:
#         validator.logger.info("Validating headers...")
#         if not validator.validate_headers():
#             validator.logger.info("Invalid headers...exiting before any further checks")
#             is_ok_to_run_validation = 0

#     if is_ok_to_run_validation:
#         validator.logger.info("Validating data...")
#         validator.validate_data()

#     # Close log handler
#     validator.logger.removeHandler(validator.handler)
#     validator.handler.close()