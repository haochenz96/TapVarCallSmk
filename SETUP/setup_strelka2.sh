if [ -f ~/.bashrc ] ; then
    . ~/.bashrc
fi

# ------> https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/quickStart.md

# change into the desired directory
cd ~/work/software

# download strelka binary
wget https://github.com/Illumina/strelka/releases/download/v2.9.2/strelka-2.9.2.centos6_x86_64.tar.bz2
# decompress
tar xvjf strelka-2.9.2.centos6_x86_64.tar.bz2 && rm strelka-2.9.2.centos6_x86_64.tar.bz2

# # run demo to check successful installation
# bash strelka-2.9.2.centos6_x86_64/bin/runStrelkaSomaticWorkflowDemo.bash
# bash strelka-2.9.2.centos6_x86_64/bin/runStrelkaGermlineWorkflowDemo.bash