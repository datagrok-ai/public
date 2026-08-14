#name: cvmBlob
#input: blob blob
#output: blob result
def cvmBlob(blob):
    return blob

#name: cvmFile
#input: file file
#output: blob result
def cvmFile(file):
    return file

#name: cvmMultiInput
#input: blob blob
#input: dataframe df
#output: string result
def cvmMultiInput(blob, df):
    return f"{len(blob)}:{df.shape[0]}:{df.shape[1]}"
