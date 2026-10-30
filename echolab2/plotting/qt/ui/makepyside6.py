#!/usr/bin/env python


import os
import stat
import sys
import logging
from PySide6 import QtCore


class makepyside6(QtCore.QObject):

    def __init__(self, make_path, pyuic):

        super(makepyside6, self).__init__()

        self.make_path = os.path.normpath(make_path)
        self.pyuic = pyuic
        
        timer = QtCore.QTimer(self)
        timer.timeout.connect(self.build_ui)
        timer.setSingleShot(True)
        timer.start(0)


    def build_ui(self):
        
         #  bump the prompt
        print()

        #  create a logger to log to the console and to a file
        logger = logging.getLogger(__name__)
        logger.propagate = False
        logger.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s : %(levelname)s : %(message)s')
        
        #  create the console logger
        consoleLogger = logging.StreamHandler(sys.stdout)
        consoleLogger.setFormatter(formatter)
        logger.addHandler(consoleLogger)
            
        logger.info("Building UI files in  " + self.make_path)

        failed = 0
        process = QtCore.QProcess()
        for name in os.listdir(self.make_path):
            source = os.path.join(self.make_path, name)
            target = None
            if source.endswith(".ui"):
                new_name = "ui_" + name.replace(".ui", ".py")
                target = os.path.join(self.make_path, new_name)

            if target is not None:
                if (not os.access(target, os.F_OK) or 
                   (os.stat(source)[stat.ST_MTIME] > os.stat(target)[stat.ST_MTIME])):
                       
                    args = ["-o", target, source]
                    process.start(self.pyuic, args)
                    
                    if not process.waitForFinished(2 * 60 * 1000):
                        logger.warning("Failed to convert " + name + ".")
                        failed += 1
                    else:
                        logger.info("converted " + source + " to " + target )
                else:
                    logger.info(new_name + " is up-to-date.")
        if failed:
            logger.warning("Errors encountered during conversion.")

        QtCore.QCoreApplication.instance().quit()
    

if __name__ == '__main__':
    
    import argparse

        #  parse the command line arguments
    parser = argparse.ArgumentParser(description='Converts one or more Qt .ui files to .py files.')
    parser.add_argument('make_path', type=str, nargs='?', default=os.getcwd())
    parser.add_argument("--pyuic_path", help="Specify the path to the pyuic script.")
    args = parser.parse_args()
    if (args.pyuic_path):
        pyuic_path = os.path.normpath(str(args.pyuic_path))
    else:
        pyuic_path = 'pyside6-uic'
    make_path = os.path.normpath(str(args.make_path))
    
    app = QtCore.QCoreApplication(sys.argv)
    form = makepyside6(make_path, pyuic_path)
    sys.exit(app.exec())
    
