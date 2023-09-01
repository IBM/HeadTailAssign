class Logger(object):
    '''A class that represents extraction of data
    Retrieved and updated from: https://stackoverflow.com/a/11325249
    
    Methods:
        write(self, obj): Write message on stdout and file.
        flush(self): Clears the internal buffer of the file.
    '''
    def __init__(self, *files):
        '''Initialize the instance of a class.

        Arguments:
        files: Files of interest.'''
        self.files = files
        for f in self.files:
            f.write("######### HEAD AND TAIL ASSIGNER LOGGER ##########")

    def write(self, obj):
        '''Write message on stdout and file.
        
        Arguments:
            obj(string): Message to be written.'''
        for f in self.files:
            f.write(obj)
            f.flush() 

    def flush(self):
        '''Clears the internal buffer of the file.'''
        for f in self.files:
            f.flush()