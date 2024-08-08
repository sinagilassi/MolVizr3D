# UTILITY CLASS
# ---------------

# import libs
import os
import json
import fnmatch
import re
import csv
import random
import numpy as np
from datetime import date
from random import randint


class Utility():

    json_csv_columns = [
        'file_name',
        'image_name_1',
        'image_name_2',
        'image_name_3',
        'mat_structure',
        'atom_numbers',
        'mat_cid',
        'mat_name',
        'mat_formula',
        'mat_smiles',
        'mat_mass',
        'atom_elements',
        'atom_atomic_number',
        'bond_numbers',
        'bond_list',
        'xyz_list',
        'xyz_center_list',
        'atom_weights',
        'mat_inchikey'
    ]

    def __init__(self):
        pass

    @staticmethod
    def CheckFileFormat(filePath):
        '''
        check file format

        args:
            filePath: file name dir

        return:
            file directory, file name, file format
        '''
        # check file exist
        if os.path.isfile(filePath):
            # file analysis
            fileDir = os.path.dirname(filePath)
            fileName = os.path.basename(filePath)
            fileFormat = os.path.splitext(filePath)[1]
            fileFormat = str(fileFormat.split(".")[-1]).lower()
            # res
            return fileDir, fileName, fileFormat
        else:
            raise Exception('file path is not valid.')

    @staticmethod
    def ReadContent(sourceContent):
        '''
        Read (load) file with respect to its format extension

        Parameters
        ----------
        sourceContent : dict
            contentFile: string content of a text file (url request)
            contentFormat: content format (sdf, json, ...)

        Returns
        -------
        file content : str
            content of a text file
        _contentFormat : str
            content format
        '''
        try:
            #
            _contentSource, _contentFormat = sourceContent.values()
            # check
            if _contentSource and _contentFormat:
                # read file
                if _contentFormat == 'sdf':
                    # string
                    fileContent = _contentSource
                elif _contentFormat == 'json':
                    # json convert to dict
                    fileContent = json.loads(_contentSource)
                elif _contentFormat == 'json-string':
                    # json convert to dict
                    fileContent = json.loads(_contentSource)

                # res
                return fileContent, _contentFormat
            else:
                raise Exception("target path is not valid.")
        except Exception as e:
            raise Exception(e)

    @staticmethod
    def OpenFile(filePath):
        '''
        Open file with respect to its format extension

        Parameters
        ----------
        filePath : str
            file path

        Returns
        -------
        file content : str
            content of a text file
        fileDir : str
            file directory
        fileName : str
            file name
        fileFormat : str
            file format
        '''
        try:
            # check
            if os.path.exists(filePath):
                # file info
                fileDir, fileName, fileFormat = Utility.CheckFileFormat(
                    filePath)

                # read a file
                with open(filePath, 'r') as f:
                    if fileFormat == 'sdf':
                        fileContent = f.read()
                    elif fileFormat == 'json':
                        fileContent = json.load(f)

                # res
                return fileContent, fileDir, fileName, fileFormat
            else:
                raise Exception("target path is not valid.")

        except Exception as e:
            raise Exception(e)

    @staticmethod
    def ListFiles(targetPath, fileExtension=''):
        '''
        list files in a target file

        Parameters
        ----------
        targetPath : str
            target path
        fileExtension : str
            file extension, default is empty

        Returns
        -------
        filesFound : list
            list of files in the target path
        '''
        try:
            # check
            if os.path.exists(targetPath):
                if not fileExtension:
                    # files
                    filesFound = os.listdir(targetPath)
                else:
                    # res
                    filesFound = []
                    for f in os.listdir(targetPath):
                        if f.endswith('.'+str(fileExtension).strip()):
                            filesFound.append(f)
                # res
                return filesFound
            else:
                raise Exception("target path is not valid.")

        except Exception as e:
            raise Exception(e)

    @staticmethod
    def SaveFile(fileContent, fileName, fileFormat, fileDir, logMessage=' file is successfully created and saved in'):
        '''
        save a file with respect to a format

        Parameters
        ----------
        fileContent : str
            file content
        fileName : str
            file name
        fileFormat : str
            file format
        fileDir : str
            file directory
        logMessage : str
            log message

        Returns
        -------
        bool
            True if the file is successfully created and saved in
        '''
        try:
            # set
            _fileFormat = str(fileFormat).lower()
            # file full name with format
            _fileFullName = fileName + f'.{_fileFormat}'
            # file full location
            fileLoc = os.path.join(fileDir, _fileFullName)

            # check
            if not os.path.isdir(fileDir):
                raise Exception('file directory is not valid.')

            # open file
            with open(fileLoc, 'w') as f:
                # check
                if _fileFormat == 'sdf':
                    # save
                    f.write(fileContent)
                    f.close()
                # FIXME
                if _fileFormat == 'json-string':
                    # save
                    json.dumps(fileContent, f, indent=5)
                    f.close()
                if _fileFormat == 'json':
                    # save ()
                    json.dump(fileContent, f, indent=5)
                    f.close()

            # log
            print(f"the {_fileFormat + logMessage} `{fileLoc}`")

            return True
        except Exception as e:
            raise Exception(e)

    @staticmethod
    def FindTargetFileInFolder(file_location, file_format):
        '''
        Find target files with respect to a specific format such as json

        Parameters
        ----------
        file_location : str
            file location
        file_format : str
            file format

        Returns
        -------
        listFiles : list
            list of target files
        '''
        try:
            # access directory
            allFiles = os.listdir(file_location)
            # file format
            fileExtension = f"*.{file_format}"
            # list
            listFiles = fnmatch.filter(allFiles, fileExtension)
            # res
            return listFiles
        except Exception as e:
            raise Exception(e)

    def SelectFile(list_file_names, file_name_ids, file_name_prefix='cid'):
        '''
        select files from a file list based on a query

        Parameters
        ----------
        list_file_names : list
            list of file names
        file_name_ids : list
            list of file ids
        file_name_prefix : str
            file name prefix

        Returns
        -------
        listFiles : list
            list of selected files
        '''
        try:
            # res
            listFiles = []
            cids = []

            # check
            if file_name_prefix == 'cid':
                # *** load only specific cids ***
                if len(file_name_ids) > 0 and len(list_file_names) > 0:
                    #
                    for i in range(len(list_file_names)):
                        _res = re.search(
                            r'(^cid_([0-9]+)(\.)([a-z0-9]*))', str(list_file_names[i]))
                        # check
                        if _res:
                            # save
                            _cid = _res.group(2)
                            cids.append(_cid)
                            # check
                            if _cid in file_name_ids:
                                _fileName = _res.group(1)
                                # save
                                listFiles.append(_fileName)
            elif file_name_prefix == 'file_name':
                # *** load only specific names ***
                if len(file_name_ids) > 0 and len(list_file_names) > 0:
                    #
                    _fileNames = [
                        item for item in file_name_ids if item in list_file_names]
                    # save
                    listFiles.extend(_fileName)
            # res
            return listFiles

        except Exception as e:
            raise Exception(e)

    @staticmethod
    def CreateCVS(csv_list, file_location, file_name='', file_format='json', save_location=''):
        '''
        Create csv file from dictionary

        Parameters
        ----------
        csv_list : list
            list of dictionary
        file_location : str
            file location
        file_name : str
            file name
        file_format : str
            file format
        save_location : str
            save location

        Returns
        -------
        _save_location : str
            save location
        csvFileName : str
            csv file name
        '''
        try:
            # check
            if len(save_location) == 0:
                _save_location = file_location
            else:
                _save_location = save_location

            # log
            _log = f"csv file is successfully created in {_save_location}"

            # csv file
            csvFileId = str(date.today()) + \
                "-" + str(random.randint(0, 256))

            # # file name
            if len(file_name) > 0:
                csvFileName = f'{file_name}.csv'
            else:
                csvFileName = f'cid_{csvFileId}.csv'

            # create path
            csvFile = os.path.join(_save_location, csvFileName)

            # write
            with open(csvFile, 'w') as f:
                writer = csv.DictWriter(
                    f, fieldnames=Utility.json_csv_columns)
                writer.writeheader()
                writer.writerows(csv_list)
                # log
                # print(_log)

            # res
            return _save_location, csvFileName
        except Exception as e:
            raise

    def ConvertCSVContentToList(csv_content):
        '''
        convert csv content (string) to list

        Parameters
        ----------
        csv_content : str
            csv string content

        Returns
        -------
        list
            list of csv content

        hint:
            list[0]: column head
            list[1]: records

        '''
        try:
            # csv column
            csv_column = []
            csv_rows = []

            #  csv list
            csv_list = csv_content.splitlines()
            # check
            if type(csv_list) is list:
                # len
                csv_list_len = len(csv_list)
                # csv column
                _column = str(csv_list[0]).replace(
                    "'", "").replace('"', '').split(',')
                csv_column.extend(_column)
                # csv rows
                _rows = str(csv_list[1]).replace(
                    "'", "").replace('"', '').split(',')
                csv_rows.extend(_rows)
            else:
                csv_list_len = -1

            # res
            return csv_list, csv_column, csv_rows, csv_list_len

        except Exception as e:
            raise Exception(e)

    @staticmethod
    def SaveNpArray(arr, file_name='', location=''):
        '''
        Save numpy array

        Parameters
        ----------
        arr : numpy array
            numpy array
        file_name : str
            file name
        location : str
            location

        '''
        try:
            # full file name
            if not file_name:
                file_name = "np_array"
            fullFileName = f"{file_name}.npy"
            # location
            fileLoc = os.path.join(location, fullFileName)

            # open
            with open(fileLoc, 'wb') as f:
                np.save(f, arr)

            print(f"save operation is done.")
        except Exception as e:
            raise Exception(e)
