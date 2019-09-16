"""
Script that checks that we have all databases and software tools
installed
"""
import os
import logging
from shutil import copyfile
from distutils.version import LooseVersion
import importlib

class Diagnostics:
    def __init__(self):
        self.logger = logging.getLogger("diagnostic_application")
        self.logger.setLevel(logging.INFO)
        formatter = logging.Formatter("%(asctime)-15s %(levelname)s %(message)s")
        file_handler = logging.FileHandler("diagnostics.log")
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)

        self.anarci = False
        self.scalop = False
    def __call__(self):
        self.logger.info("\tWriting DIAGNOSTICS log\n")    
       
        self.check_anarci()
        self.check_scalop()
        self.check_fread()
        self.check_fread_db()
        self.check_dependencies()
        print "Diagnostics log is found in diagnostics.log"

    def check_anarci(self):

        # Importing anarci
        try:
            from anarci import run_anarci
            self.logger.info(" Successfully imported: anarci")
        except ImportError:
            self.logger.error(" Cannot import: anarci")

        # Importing anarci germlines
        try:
            from anarci.germlines import all_germlines
            self.logger.info(" Successfully imported: anarci germlines")
            self.anarci = True
        except ImportError:
            self.logger.error(" Cannot import: anarci germlines")
        # Imporint anarci Accept class
        try:
            from saab_plus.aboss_utils.region_definitions import Accept
            self.logger.info(" Successfully imported: anarci Accept")
        except ImportError:
            self.logger.error(" Cannot import: anarci Accept")

    def check_scalop(self):
        
        try:
            import scalop
            self.logger.info(" Successfully imported: scalop")
            self.scalop = True
        except:
            self.logger.error(" Cannot import: scalop")

        if self.scalop:
            try:
                from scalop.inhouse_predict import _assign
                self.logger.info(" Successfully imported: scalop assing")
            except ImportError:
                # copy anarci germline to scalop 
                self.copy_germlines()

                try:
                    from scalop.inhouse_predict import _assign
                    self.logger.info(" Successfully imported: scalop assing")
                except:
                    self.logger.error(" Cannot import: scalop")

    def copy_germlines(self):
        """
        Copying anarci germlines to scalop directory, since
        scalop does not come with pre-installed germlines
        """
        self.logger.warning(" scalop does not have germline!!!")
        import anarci
        anarci_loc = os.path.dirname(anarci.__file__)
        anarci_germ = os.path.join(anarci_loc, "germlines.py")
        if not os.path.isfile(anarci_germ):
            self.logging.error(" Cannot locate: anarci germlines.py")
        
        import scalop
        scalop_loc = os.path.dirname(scalop.__file__)
        scalop_germ =  os.path.join(scalop_loc, "anarci/germlines.py")
        copyfile(anarci_germ, scalop_germ)
        self.logger.info(" anarci germlines were copied to scalop")

    def check_fread(self):
       
        try:
            import code.Alignment.FREAD.pyfread_api as fread
            self.logger.info(" Successfully imported: FREAD")
            esst = fread.esst_path

            if os.path.isfile(esst):
                self.logger.info(" Successfully imported: FREAD ESS table")
            else:
                self.logger.error(" FREAD ESS table is not found")
        except ImportError:
            self.logger.error(" Cannot import: FREAD")

        try:
            import code.Alignment.FREAD.prosci
            self.logger.info(" Successfully imported: prosci module")
        except ImportError:
            self.logger.error(" Cannot import: prosci module")

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
            self.logger.info(" Successfully imported: Common module")
        except ImportError:
            self.logger.error(" Cannot import: Common module")
        try:
            from code.Common.Common import pdb_length_location, template_db,\
                                            fread_db, numbered_datasets_location,\
                                            clusters, numbered_germlines
            # checking pdb_length
            if not self.file_exist( pdb_length_location ):
                self.logger.error(" PDB template info file is not found!")
            else:
                self.logger.info("         PDB template info file is located: OK")
            
            if not self.file_exist(clusters):
                self.logger.error(" CDR-H3 cluster mapping file is missing!")
            else:
                self.logger.info("    CDR-H3 cluster mapping file is located: OK")

            # locating PDB frameworks
            if not self.dir_exist(template_db):
                self.logger.error(" Directory with PDB frameworks not found")
            else:
                if not os.listdir(template_db):
                    self.logger.error(" Directory with PDB frameworks is empty")
                else:
                    self.logger.info("  Directory with PDB frameworks is located: OK")

            # locating fread pdb templates
            if not self.dir_exist(fread_db):
                self.logger.error( "FREAD PDB template directory not found")
            else:
                self.logger.info(" Directory with FREAD templates is located: OK")
                self.logger.info(" Number of FREAD templates found: {0}".format(self.count_templates(fread_db, ".atm.gz")))

            if not self.dir_exist( os.path.join( numbered_datasets_location,'sabdab', "IMGT")):
                self.logger.error(" Directory with numbered PDB frameworks is not found")
            else:
                self.logger.info(" Directory with numbered PDB frameworks: OK")
                if not os.listdir( os.path.join( numbered_datasets_location,'sabdab', "IMGT") ):
                    self.logger.error(" Directory with numbered PDB frameworks: empty!!!")
        except ImportError:
            self.logger.error(" Cannot import PATHs from Common module")

    def check_dependencies(self):
        """
        Checking if a user has Biopython and pandas
        """
        dependencies = [("pandas", "0.24.2"),
                        ("Bio", "1.70")]
        for lib, version in dependencies:
            try:
#                globals()[lib] = importlib.import_module(lib)
                 imported_lib = importlib.import_module(lib)
                 imported_version = imported_lib.__version__ 
                 if LooseVersion(imported_version) < LooseVersion(version):
                    self.logger.warning(" Package version is older than recommended for : {0}".format(imported_lib.__name__)) 
                    self.logger.warning(" Imported version is {0}, recomimended {1}".format( imported_version, version))
                    print "WARGNING: {0} version is older than recommended".format(imported_lib.__name__)
                    
                 else:
                    self.logger.info(" Imported {0}, Version {1}".format( imported_lib.__name__, imported_lib.__version__))

            except:
                self.logger.error(" Cannot import package: ", lib)
                print "Cannot import: ", lib

if __name__ == "__main__":
    diagnostics = Diagnostics()
    diagnostics()
