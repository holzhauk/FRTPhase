class SPhaseFileError(Exception):
    pass

class SPhaseFileIOError(SPhaseFileError):
    def __init__(self, file_name, message):
        self.file_namee = file_name
        self.message = message

class SPhaseFileWrongFormat(SPhaseFileError):
    def __init__(self, file_name, message):
        self.file_name = file_name
        self.message = message

class SPhaseFileVersionConflict(SPhaseFileError):
    def __init__(self, file_name, message):
        self.file_name = file_name
        self.message = message

class SPhaseFileModelConflict(SPhaseFileError):
    def __init__(self, file_name, message):
        self.file_name = file_name
        self.message = message

class SPhaseFileClassConflict(SPhaseFileError):
    def __init__(self, file_name, message):
        self.file_name = file_name
        self.message = message

class SPhaseFileDimError(SPhaseFileError):
    def __init__(self, file_name, message):
        self.file_name = file_name
        self.message = message