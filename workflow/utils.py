def parse_for_shell_args(input_dict):
    '''
    General function to parse a python dictionary into bash arguments
    '''
    args = ""
    for key, value in input_dict.items():
        args += f"-{key} {value} "

    return args

def parse_bcf_filters(filter_list):
    '''
    Parse a Python list into bcftools filter arguments
    
    '''
    filter_string = " & ".join(filter_list)
    return filter_string

