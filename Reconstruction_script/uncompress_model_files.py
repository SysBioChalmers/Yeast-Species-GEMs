# -*- coding: utf-8 -*-
# -*- python 3 -*-
# -*- hongzhong Lu -*-

# parse the rhea database

import os    ##for directory
# uncompress
os.system("gunzip -k xxxx/GitHub/Yeast-Species-GEMs/Reconstruction_script/ModelFiles/*.gz")

# compress
os.system("cd GitHub/Yeast-Species-GEMs/Reconstruction_script/ModelFiles/")
os.system("tar -czvf xml.tar.gz xml")
os.system("tar -czvf txt.tar.gz txt")

