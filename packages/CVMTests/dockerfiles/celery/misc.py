import time

#name: cvmApiKey
#output: string result
def cvmApiKey(**kwargs):
    return str(kwargs.get('USER_API_KEY', ''))

#name: cvmCancel
#input: int seconds
#output: string result
def cvmCancel(seconds):
    time.sleep(seconds)
    return 'done'

#name: cvmTwoOutputs
#output: int a
#output: int b
def cvmTwoOutputs():
    return 1, 2
