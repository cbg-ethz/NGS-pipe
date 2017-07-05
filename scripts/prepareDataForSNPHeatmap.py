import sys

inFileName = sys.argv[1]
outFileName = sys.argv[2]


outFile = open(outFileName, 'w')
inFile = open(inFileName, 'r')

for line in inFile:
    if line.startswith("##"):
        continue

    lineSplit = line.strip().split("\t")
    if line.startswith("#"):
        for i in range(9, len(lineSplit)):
            outFile.write(lineSplit[i] + "\t")
        outFile.write("\n")
        continue

    formatSplit = lineSplit[8].split(":")
    try:
        gtIndex = formatSplit.index("GT")
    except:
        print("No GT field for position: " + lineSplit[0] + "\t" + lineSplit[1])
        continue

    for i in range(9, len(lineSplit)):
        sampleSplit = lineSplit[i].split(":")
        gt = sampleSplit[gtIndex]
        if gt == "./.":
            gt = "NA"
        outFile.write(gt + "\t")
    outFile.write("\n")


