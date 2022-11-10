# (c) 2019. All Rights Reserved.
# Code written by: 
#       Maxim
#       Hugo
#       Sean Powell (sean.powell@imba.oeaw.ac.at)

def read_samples_file(samples_file, include_column = None):
    (_, file_ext) = os.path.splitext(samples_file)

    if file_ext == '.csv':
        df = pd.read_csv(samples_file)
    elif file_ext == '.xls':
        df = pd.read_excel(samples_file)
    elif file_ext == '.xlsx':
        df = pd.read_excel(samples_file)
    else:
        return None

    if include_column is None:
        return df

    return df[df[include_column] == True]


OUTPUT_EXT = '.{resolution}{output_ext}'

def get_cooler_files(grouped_dict, resolution, output_ext):
    cooler_files = list()
    for key in grouped_dict.keys():
        if isinstance(key, str) or isinstance(key, int):
            cooler_file_name = '{}{}'.format(key, OUTPUT_EXT.format(resolution = resolution, output_ext = cooler_ext)) 
        else:
            arr = [str(k) for k in key if not k == '']       
            cooler_file_name = '{}{}'.format('_'.join(arr), OUTPUT_EXT.format(resolution = resolution, output_ext = cooler_ext))
            
        cooler_file = os.path.join(cooler_dir, cooler_file_name)
        if not os.path.exists(cooler_file):
            sys.stderr.write('warning: merged cooler missing {}\n'.format(cooler_file))
            continue

        cooler_files.append(cooler_file)

    return cooler_files

def group_by_tags(df, tag_names):
    df = df.fillna('')
    df = df.sort_values(tag_names)
    grouped_dict = df.groupby(tag_names).groups

    return grouped_dict

def get_chunks(l, chunks):
    a = [t for t in zip(*[iter(l)] * chunks)]: # split the array l in to chunk sized portions
        
    if not len(l) % chunks == 0: # add the final chunk
        a.append(tuple(l[-(len(l) % chunks):]))

    return a