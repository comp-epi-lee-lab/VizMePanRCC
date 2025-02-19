import pandas as pd

def read_coords_file(file_path):
    df = pd.read_csv(file_path, sep='\t', index_col='Composite_Element_REF')
    return df

def get_cg_site_details(df, cg_sites):
    gene_symbols = []
    chromosomes = []
    genomic_coordinates = []

    for cg_site in cg_sites:
        if cg_site in df.index:
            row = df.loc[cg_site]
            gene_symbols.append('' if pd.isna(row['Gene_Symbol']) else row['Gene_Symbol'])
            chromosomes.append('' if pd.isna(row['Chromosome']) else row['Chromosome'])
            genomic_coordinates.append('' if pd.isna(row['Genomic_Coordinate']) else row['Genomic_Coordinate'])
        else:
            gene_symbols.append(None)
            chromosomes.append(None)
            genomic_coordinates.append(None)

    return gene_symbols, chromosomes, genomic_coordinates

def main():
    file_path = r'data/KIPAN.hm450.coords' 
    df_coords = read_coords_file(file_path)

    cg_sites = [
    'cg04096619',
    'cg05343811'
    ]

    gene_symbols, chromosomes, genomic_coordinates = get_cg_site_details(df_coords, cg_sites)

    print("Gene Symbols:")
    for symbol in gene_symbols:
        print(symbol)

    # print("\nChromosomes:")
    # for chromosome in chromosomes:
    #     print(chromosome)

    # print("\nGenomic Coordinates:")
    # for coordinate in genomic_coordinates:
    #     print(coordinate)


if __name__ == '__main__':
    main()