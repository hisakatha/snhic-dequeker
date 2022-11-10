
import os
import glob
import click
import pandas as pd

BAM_FOLDER = '/groups/tachibana/sean/projects/johanna/snhic/mouse/dmog.pairtools/data/bam'
CSV_FILE = '/groups/tachibana/sean/projects/johanna/snhic/mouse/dmog.pairtools/johanna.dmog.maternal.paternal.csv'
BAM_FILE = '*{}*.bam'
BAM_FILES = '*.bam'
def main(bam_folder, csv_file):
    df = pd.read_csv(csv_file)
    df = df[df['include'] == True]

    sampleIDs = df.cellID.values

    sampleIDs.sort()
    arrBamFiles = list()
    for sampleID in sampleIDs:
        bam_files = glob.glob(os.path.join(bam_folder, BAM_FILE.format(sampleID)))
        if len(bam_files) == 0:
            print('missing: cellID {}'.format(sampleID))

        elif len(bam_files) > 1:
            print('multiple: cellID {}'.format(sampleID))
            print(bam_files)

        else:
            arrBamFiles.append(bam_files[0])

    setExperimentFiles = set(arrBamFiles)
    setBamsFilesListed = set(glob.glob(os.path.join(bam_folder, BAM_FILES)))

    setBamFilesExtra = setBamsFilesListed.difference(setExperimentFiles)
    bamFilesExtra = list(setBamFilesExtra)
    bamFilesExtra.sort()
    for bam_file in bamFilesExtra:
        print(bam_file)
        


if __name__ == '__main__':
    main(BAM_FOLDER, CSV_FILE)