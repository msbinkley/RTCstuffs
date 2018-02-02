# RTCstuffs

# Description
    - A git repository to hold the code to run RTC stuffs. 
    - This will allow local development and testing, which we can then run in full force on a server. 
    - We will make sure to that all of the variables referenced will be easily changed for the server (or any other location where the analysis can be run)

# How to run. 
    1. Edit your paramfile to match your directory structure.
        - binkley_local for Binkley's MacBook Pro
    2. Use the sample jupyter notebook for how to load the param_file data. 

    3. Keep your jupyter notebooks in git directory so we can share our work and progress. 

    4. The variables in the paramfiles can be accessed by: 
        paramDict["   "]

# PARAMETER FILE VARIABLES
    - tmpDir: A directory where the temporary files are kept. 
    - gtexDir: A directory for the gtex data
    - vcfDir: A directory for the gtex vcf directory (should be in gtexDir)

# Notes: 
    - Keep files in the git directory small.  (<10mb)
        - Git uploads can take a long time when there are large files. 
        - Keep the large files (like VCF or GTEx data in a separate directory outside of the git directory. 
        - The locations of the large files can be referenced in the param files. 
