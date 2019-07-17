"""
Script that checks that we have all databases and software tools
install on the system
"""
import logging
import os

class Diagnostics:

    def __init__(self):
        pass    

    def check_anarci(self):

        # Importing anarci
        try:
            from anarci import run_anarci
            print "Successfully imported: anarci"
        except ImportError:
            print "Cannot import: anarci"

        # Importing anarci germlines
        try:
            from anarci.germlines import all_germlines
            print "Successfully imported: anarci germlines"
        except ImportError:
            print "Cannot import: anarci germlines"

    def check_scalop(self):
        
        try:
            from scalop.inhouse_predict import _assign
            print "Successfully imported: scalop"
        except ImportError:
            print "Cannot import: scalop"

    def check_fread(self):
       
        try:
            import code.Alignment.FREAD.pyfread_api as fread
            print "Successfully imported: FREAD"
            where_i_am = fread.where_am_i
            esst = os.path.join(where_i_am, "esst.txt")
            if os.path.isfile(esst):
                print "Successfully imported: FREAD ESS table"
            else:
                print "FREAD ESS table is not found"
        except ImportError:
            print "Cannot import: FREAD"

        try:
            import code.Alignment.FREAD.prosci
            print "Successfully imported: prosci module"
        except ImportError:
            print "Cannot import: prosci module"

    def file_exist(self, file_name):
        if os.path.isfile(file_name):
            return True
        return False
    
    def dir_exist(self, dir_name):
        if os.path.isdir(dir_name):
            return True
        return False

    def count_templates(self, dir_name, extension):
        total = 0
        for _, _, pdbs in os.walk(dir_name):
            for pdb in pdbs:
                if extension in pdb:
                    total +=1
        return total
    
    def check_fread_db(self):

        try:
            import code.Common.Common
            print "Successfully imported: Common module"

        except ImportError:
            print "Cannot import: Common module"

        try:
            from code.Common.Common import pdb_length_location, template_db,\
                                            fread_db, numbered_datasets_location
            # checking pdb_length
            if not self.file_exist( pdb_length_location ):
                print "PDB template info file is not found!"
            else:
                print "           PDB template info file is found: {0}".format(pdb_length_location)
            
            # locating PDB frameworks
            if not self.dir_exist(template_db):
                print "Directory with PDB frameworks not found"
            else:
                if not os.listdir(template_db):
                    print "Directory with PDB frameworks is empty"
                else:
                    print "  Directory with PDB frameworks is located: {0}".format(template_db)

            # locating fread pdb templates
            if not self.dir_exist(fread_db):
                print "FREAD PDB template directory not found"
            else:
                print " Directory with FREAD templates is located: {0}".format(fread_db)
                print " Number of FREAD templates found: {0}".format(self.count_templates(fread_db, ".atm.gz"))

            if not self.dir_exist( os.path.join( numbered_datasets_location,'sabdab', "IMGT")):
                print "Directory with numbered PDB templates is not found"
            else:
                print "     Directory with numbered PDB templates: {0}".format(os.path.join( numbered_datasets_location,'sabdab', "IMGT"))
                if not os.listdir( os.path.join( numbered_datasets_location,'sabdab', "IMGT") ):
                    print "Directory with numbered PDB templates: empty!!!"
            
        except ImportError:
            print "Cannot import PATHs from Common module"

if __name__ == "__main__":
    
    diagnostics = Diagnostics()
    diagnostics.check_anarci()
    diagnostics.check_scalop()
    diagnostics.check_fread()
    diagnostics.check_fread_db()
