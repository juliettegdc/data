#!/usr/bin/env python
#
#    Copyright (C) 2013 Alexandros Avdis and others. See the AUTHORS.md file for a full list of copyright holders.
#
#    This file is part of QMesh.
#
#    QMesh is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    QMesh is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with QMesh.  If not, see <http://www.gnu.org/licenses/>.

from setuptools.command.install import install as _install

def main():
    from setuptools import setup
    import os
    import subprocess
    # Put files containing git sha key, version, license, authors list and
    # README in the right place for the egg-file.
    if not os.path.isfile('qmesh/GIT_SHA_KEY'):
        git_sha_key_copied = subprocess.call(['cp','GIT_SHA_KEY','qmesh/'])
        if git_sha_key_copied != 0:
            try:
                git_sha_key = subprocess.check_output(['git','rev-parse','HEAD'])
            except:
                git_sha_key = 'Could not obtain git sha key.'
            git_sha_key_file = open('qmesh/GIT_SHA_KEY','w')
            git_sha_key_file.write(git_sha_key.strip())
            git_sha_key_file.close()
    license_copied = subprocess.call(['cp','LICENSE','qmesh/'])
    authors_copied = subprocess.call(['cp','AUTHORS.md','qmesh/'])
    readme_copied = subprocess.call(['cp','README.md','qmesh/'])
    version_copied = subprocess.call(['cp','VERSION','qmesh/'])
    #Read version from file
    version_file = open('qmesh/VERSION','r')
    qmesh_version_string = version_file.readline().strip()
    version_file.close()

    try:
      destdir = os.environ["DESTDIR"]
    except KeyError:
      destdir = "/usr/share/"
    try:
        set
    except NameError:
        from sets import Set
  
    setup(
          name='qmesh',
          version=qmesh_version_string,
          description = "Finite Element meshes from GIS data.",
          author = "The QMesh Development Team.",
          author_email = "develop@qmesh.org",
          url = "http://www.qmesh.org",
          download_url = 'https://bitbucket.org/qmesh-developers/qmesh/commits/tag/1.0.0',
          packages = ['qmesh',
                      'qmesh.lib',
                      'qmesh.vector',
                      'qmesh.mesh',
                      'qmesh.raster',
                      'qmesh.publish',
                     ],
          package_dir = {
              'qmesh': 'qmesh',
              'qmesh.lib':'qmesh/lib',
              'qmesh.vector':'qmesh/vector',
              'qmesh.meshg':'qmesh/mesh',
              'qmesh.raster':'qmesh/raster',
              'qmesh.publish':'qmesh/publish',
              },
          provides=['qmesh'],
          install_requires=['GDAL', 'sword2', 'restkit', 'GitPython', 'requests'],
          package_data = {'qmesh':['VERSION','GIT_SHA_KEY','LICENSE','AUTHORS.md','README.md']},
          license='GPLv3',
          test_suite = "tests",
          keywords = ['GIS', 'mesh generation'],
          cmdclass={'install': install}
        )

    # A little clean-up
    #subprocess.call(['rm','qmesh/VERSION','qmesh/GIT_SHA_KEY','qmesh/LICENSE','qmesh/AUTHORS.md','qmesh/README.md'])

def _check_gmsh():
    '''Check gmsh is installed, by trying to find out installed version'''
    import logging
    import subprocess
    logging.basicConfig(format='%(message)s',level=logging.DEBUG)
    this_pid = subprocess.os.getpid()
    try:
        gmsh_stdoutFileName = '/tmp/gmsh_stdout_'+str(this_pid)
        gmsh_stdout = open(gmsh_stdoutFileName,'w')
        gmsh_found = subprocess.Popen(['gmsh','--version'],stdout = gmsh_stdout, stderr = gmsh_stdout)
        gmsh_found.wait()
        gmsh_stdout.close()
        if gmsh_found.returncode != 0:
            msg = 'Could not find a gmsh installation.'
            raise Exception(msg)
        else:
            gmsh_file = open(gmsh_stdoutFileName,'r')
            line = gmsh_file.next()
            logging.info("Found gmsh version "+line.strip())
    except:
        logging.error(msg)
        raise

def _check_qgis():
    '''Check gmsh is installed, by trying to find out installed version and reading in a shapefile'''
    import logging
    logging.basicConfig(format='%(message)s',level=logging.DEBUG)
    try:
        import qgis.core
        qgs = qgis.core.QgsApplication([], False)
        qgis_install_path='/usr'
        qgs.setPrefixPath(qgis_install_path, True)
        qgs.initQgis()
        logging.info('Found QGIS version '+qgis.core.QGis.QGIS_VERSION)
    except:
        raise Exception('Could not find a qgis installation at /usr.')

def _check_pyrdm():
    '''Check if PyRDM is installed. Fetch and install if not.'''
    import logging
    import subprocess
    logging.basicConfig(format='%(message)s',level=logging.DEBUG)
    try:
        import pyrdm
        logging.info('Found PyRDM (at'+pyrdm.__path__[0]+').')
    except:
        logging.warning('Could not find PyRDM. Attempting to fetch from remote repository and install...')
        try: #Try downloading
            pyrdm_download = subprocess.Popen(['wget','https://github.com/pyrdm/pyrdm/archive/master.tar.gz'])
            pyrdm_download.wait()
            msg = 'Could not fetch PyRDM repository.'
            if pyrdm_download.returncode != 0:
                raise Exception(msg)
        except:
            logging.error(msg, exc_info = True)
            raise
        try: #Try installing
            pyrdm_expand = subprocess.Popen(['tar','-xzf','master.tar.gz'])
            pyrdm_expand.wait()
            if pyrdm_expand.returncode == 0:
                pyrdm_install = subprocess.Popen(['make','install'], cwd='pyrdm-master')
                pyrdm_install.wait()
                if pyrdm_install.returncode == 0:
                    pyrdm_clean = subprocess.Popen(['rm','-rf','master.tar.gz','pyrdm-master'])
                    pyrdm_clean.wait()
                    if pyrdm_clean.returncode == 0:
                        return
            msg = 'Could not install PyRDM.'
            raise Exception(msg)
        except:
            logging.error(msg, exc_info=True)
            raise
        raise
            

def _check_GFD_basisChangeTools():
    '''Check if GFD_basisChangeTools is installed. Fetch and install if not.'''
    import logging
    import subprocess
    logging.basicConfig(format='%(message)s',level=logging.DEBUG)
    try:
        import GFD_basisChangeTools
        logging.info('Found GFD_basisChangeTools (at '+GFD_basisChangeTools.__path__[0]+').')
    except:
        logging.warning('Could not find GFD_basisChangeTools. Attempting to fetch from remote repository and install...')
        try: #Try downloading
            GFD_bCT_download = subprocess.Popen(['wget','https://github.com/AlexandrosAvdis/GFD_basisChangeTools/archive/master.tar.gz'])
            GFD_bCT_download.wait()
            if GFD_bCT_download.returncode != 0:
                msg = 'Could not fetch GFD_basisChangeTools repository.'
                raise Exception(msg)
        except:
            logging.error(msg, exc_info = True)
            raise
        try: #Try installing
            GFD_bCT_expand = subprocess.Popen(['tar','-xzf','master.tar.gz'])
            GFD_bCT_expand.wait()
            if GFD_bCT_expand.returncode == 0:
                GFD_bCT_install = subprocess.Popen(['python','setup.py','install'],cwd='GFD_basisChangeTools-master')
                GFD_bCT_install.wait()
                if GFD_bCT_install.returncode == 0:
                    GFD_bCT_clean = subprocess.Popen(['rm','-rf','master.tar.gz','GFD_basisChangeTools-master'])
                    GFD_bCT_clean.wait()
                    if GFD_bCT_clean.returncode == 0:
                        return
            msg = 'Could not install GFD_basisChangeTools.'
            raise Exception(msg)
        except:
            logging.error(msg, exc_info=True)
            raise
        raise

class install(_install):
    def run(self):
        self.execute(_check_gmsh, (), msg='Looking for gmsh...')
        self.execute(_check_qgis, (), msg='Looking for QGIS...')
        self.execute(_check_pyrdm, (), msg='Looking for PyRDM...')
        self.execute(_check_GFD_basisChangeTools, (), msg='Looking for GFD_basisChangeTools...')
        _install.run(self)

if __name__=='__main__':
    main()
