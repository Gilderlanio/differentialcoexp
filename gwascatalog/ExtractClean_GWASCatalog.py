def read_rows(filename):
    print("Reading GWAS Catalog rows...")
    rows = []
    try:
        with open(filename, 'r') as f:
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
        print("Read GWAS Catalog rows sucessfuly!")

def write_rows_to_file(rows, output_path):
    try:
        with open(output_path, 'w') as o:
            for row in rows:
                row = [str(elem) for elem in row]
                o.write('\t'.join(row))
                o.write('\n')
            o.close()
    except IOError:
        print('Can\'t write data on file.')


def preprocess_gwascatalog(gwascatalog_rows):

    gwas_preprocessed = set()
    snps = set()

    for row in gwascatalog_rows:

        pmid = row[1]
        author = row[2] + ", " + row[3]
        phenotype = str(row[7]).lower()
        region = row[10]
        chr = row[11]
        chr_pos = row[12]
        reported_gene = row[14].upper()
        snp = row[20]
        pvalue = row[27]
        odds = row[30]

        if phenotype != "" and snp != "" and snp.count("rs") == 1:
            if snp.count("-") == 1:
                snp_info = snp.split('-')
                rs = snp_info[0]
                if rs[-1] == " ":
                    rs = snp[0:-1]
                allele = str(snp_info[1]).replace(" ", "")
            else:
                rs = snp_info[0]
                if rs[-1] == " ":
                    rs = rs[0:-1]
                allele = "?"
            if region == "":
                region = "?"
            if chr == "":
                chr = "?"
            t = (pmid, author, phenotype, chr, chr_pos, region, reported_gene, rs, allele.upper(), pvalue, odds)
            gwas_preprocessed.add(t)
            snps.add((chr, chr_pos, rs, allele.upper()))
    snps = list(snps)
    snps.insert(0, ["CHROM", "POS", "SNP"])
    return gwas_preprocessed, snps

if __name__ == '__main__':

    import sys
    file = sys.argv[1]
    gwas_rows =  read_rows(file)
    gwas_rows, snps  = preprocess_gwascatalog(gwas_rows)
    write_rows_to_file(gwas_rows, file + "_PreProcessado.txt")
    write_rows_to_file(snps, file + "_SNPs.txt")