# RTCstuffs

# Description
    - A git repository to hold the code to run RTC stuffs. 
    - This will allow local development and testing, which we can then run in full force on a server. 
    - We will make sure to that all of the variables referenced will be easily changed for the server (or any other location where the analysis can be run)

# How to run. 
    1. Create your own paramfile to match your directory structure.
        - binkley_local for Binkley's MacBook Pro
    2. Load the paramfile in the top of the jupyter notebook by:

    3. Access the functions by appending the following to the top of your jupyter notebook as well:

    4. Keep your jupyter notebooks in git directory so we can share our work. 

    5. The variables in the paramfiles can be accessed by: 
        paramDict["   "]

# PARAMETER FILE VARIABLES
    - largeTmpDir: A directory where the temporary files are kept. 

# Notes: 
    - Keep files in the git directory small.  (<10mb)
        - Git uploads can take a long time when there are large files. 
        - Keep the large files (like VCF or GTEx data in a separate directory outside of the git directory. 
        - The locations of the large files can be referenced in the param files. 
