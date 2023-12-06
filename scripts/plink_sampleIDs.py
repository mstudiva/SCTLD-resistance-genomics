with open('ofav_passing.fam', 'r') as fam_file:
    for line in fam_file:
        sample_id = line.split()[1]  # Assuming the second column is the sample ID
        print(sample_id)
