__author__ = 'GilAraujo'

def read_rows_from_tab_file(filename, header):
    """
    Read rows from tab delimited file.
    Returns an array/list with rows.

    """
    rows = []
    try:
        with open(filename, 'rU') as f:
            if header:
                next(f)
            for line in f:
                line = line.strip('\n')
                line = line.strip('\r')
                line = line.split('\t')
                rows.append(line)
            f.close()
        return rows
    except IOError:
        print("Error: can\'t read file.")
    else:
        print("Read rows sucessfuly!")

def read_row_pandas(filename, separator= "\t"):
    """
    Read rows from tab delimited file.
    Returns a pandas dataframe.
    """
    import pandas
    try:
        df = pandas.read_csv(filename, separator, low_memory=False)
        return df
    except IOError:
        raise Exception("Error: can\'t read file.")
    else:
        print("Read rows sucessfuly!")
