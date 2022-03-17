def make_all_input_mutect():
    # output is 1 vcf per patient, 
    patients = set()
    with open(config["samples"], "r") as f:
        for l in f:
            if not l.startswith("PATIENT_ID"):
                info = l.strip().split() 
                patients.add( info[0] )
        #paths =  [  os.path.join( "mutect", p, f"unfiltered_{i}.vcf") for p in patients for i in range(1,24) ]
        #paths =  [  os.path.join( "mutect", p, "read-orientation-model.tar.gz") for p in patients ]
        #paths =  [  os.path.join( "mutect", p, "filtered.vcf") for p in patients ]
        paths =  [  os.path.join( "mutect", p, "filtered_PASS_biallelic.vcf.gz") for p in patients ]
        print(paths)
        return paths

def make_all_input_intersect():
    # output is 1 vcf per patient, 
    patients = set()
    with open(config["samples"], "r") as f:
        for l in f:
            if not l.startswith("PATIENT_ID"):
                info = l.strip().split() 
                patients.add( info[0] )
        #paths =  [  os.path.join( "mutect", p, f"unfiltered_{i}.vcf") for p in patients for i in range(1,24) ]
        #paths =  [  os.path.join( "mutect", p, "read-orientation-model.tar.gz") for p in patients ]
        #paths =  [  os.path.join( "mutect", p, "filtered.vcf") for p in patients ]
        paths =  [  os.path.join(f"intersect/{p}/intersection.vcf.gz") for p in patients ]
        print(paths)
        return paths