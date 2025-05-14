import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import sys

def load_degs(file_path):
    """Load DEG gene IDs from Dataset S8.xlsx (2 columns: YPK, OPK)"""
    df = pd.read_excel(file_path, engine='openpyxl')
    ypk_genes = set(df["YPK"].dropna().astype(str))
    opk_genes = set(df["OPK"].dropna().astype(str))
    return ypk_genes, opk_genes

def load_region_genes(file_path, region):
    """Load gene IDs for the specified region (Promoter, Exon, Intron)"""
    df = pd.read_excel(file_path, engine='openpyxl')
    return set(df[region].dropna().astype(str))

def plot_venn3_custom(ypk_set, region_set, opk_set, region_label, output_path):
    """Plot Venn diagram with left=YPK, center=Region, right=OPK in correct color order"""
    plt.figure(figsize=(6, 6))
    venn = venn3([ypk_set, region_set, opk_set],
                 set_labels=("YPK", region_label, "OPK"))

    # Color patches correctly: '100' = YPK only, '010' = Region only, '001' = OPK only
    if venn.get_patch_by_id('100'):
        venn.get_patch_by_id('100').set_color('#b5cf49')  # YPK only
        venn.get_label_by_id('100').set_color('#b5cf49')
    if venn.get_patch_by_id('010'):
        venn.get_patch_by_id('010').set_color('gray')     # Region only
        venn.get_label_by_id('010').set_color('gray')
    if venn.get_patch_by_id('001'):
        venn.get_patch_by_id('001').set_color('#c38452')  # OPK only
        venn.get_label_by_id('001').set_color('#c38452')

    plt.title(f"Venn Diagram: YPK, {region_label}, OPK")
    plt.savefig(output_path)
    plt.close()
    print(f"Saved: {output_path}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python venn_regions_s8.py <Dataset_S8.xlsx> <Dataset_S7.xlsx> <output_prefix>")
        sys.exit(1)

    deg_file = sys.argv[1]
    region_file = sys.argv[2]
    output_prefix = sys.argv[3]

    # Load DEG gene sets
    ypk_genes, opk_genes = load_degs(deg_file)

    # Iterate through regions and draw customized Venn diagrams
    for region in ["Promoter", "Exon", "Intron"]:
        region_genes = load_region_genes(region_file, region)
        output_path = f"{output_prefix}_{region.lower()}.pdf"
        plot_venn3_custom(ypk_genes, region_genes, opk_genes, region, output_path)