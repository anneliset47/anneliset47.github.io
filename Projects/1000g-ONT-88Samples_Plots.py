import json
import os
import pysam  # Used to read VCF files
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd

# File paths
VCF_PATH = "/Users/annelisethorn/Documents/Anschutz/Datasets/1000g-ONT-88Samples.vcf.gz"
JSON_PATH = "/Users/annelisethorn/Documents/Anschutz/Datasets/STRchive-loci.json"
OUTPUT_DIR = "/Users/annelisethorn/Documents/Anschutz/Plots/1000g-ONT-88Samples_tr_plots"

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Load JSON metadata for loci
with open(JSON_PATH, "r") as file:
    loci_data = json.load(file)

# Function to annotate the plot with range lines and labels
def add_range_marker_or_line(fig, x0, x1, label_id, label_positions,
                              chart_width_px, x_span,
                              base_y_line=1.02, base_y_label=1.1,
                              bump_y=0.08, bump_px_threshold=60):
    
    # Adds either a horizontal line or vertical line (if a single point) to the plot
    # with a label above it. Bumps labels if overlapping.
    
    x_center = (x0 + x1) / 2
    pixels_per_unit = chart_width_px / x_span

    # Bump logic based on visual overlap in pixel space
    level = 0
    for other_x in label_positions.values():
        if abs(other_x - x_center) * pixels_per_unit < bump_px_threshold:
            level += 1

    label_positions[label_id] = x_center
    y_line = base_y_line + bump_y * level
    y_label = base_y_label + bump_y * level

    label = f"{label_id.capitalize()}" if x0 == x1 else f"{label_id.capitalize()}"

    # Add label above the line
    fig.add_annotation(
        x=x_center, y=y_label,
        text=label, showarrow=False,
        xref="x", yref="paper",
        font=dict(color="black", size=12)
    )

    # Add line (dashed if a single point)
    if x0 == x1:
        fig.add_shape(
            type="line", x0=x0, x1=x0, y0=0, y1=1,
            xref="x", yref="paper",
            line=dict(color="black", width=2, dash="dot")
        )
    else:
        fig.add_shape(
            type="line", x0=x0, x1=x1, y0=y_line, y1=y_line,
            xref="x", yref="paper",
            line=dict(color="black", width=2)
        )

# Read the VCF
vcf_in = pysam.VariantFile(VCF_PATH)

for record in vcf_in.fetch():
    chrom = record.chrom
    pos = record.pos

    # Find the matching locus in the JSON file
    matched_locus = next(
        (entry for entry in loci_data
         if chrom == entry["chrom"] and entry["start_hg38"] <= pos <= entry["stop_hg38"]),
        None
    )
    if not matched_locus:
        continue

    motif = matched_locus.get("reference_motif_reference_orientation")
    motif_length = len(motif) if motif else 1

    # Extract allele lengths and convert to repeat counts
    repeat_counts = []
    for sample in record.samples.values():
        al_lengths = sample.get("AL")
        if al_lengths:
            repeat_counts.extend([int(length) // motif_length for length in al_lengths if length is not None])

    if not repeat_counts:
        continue

    # Build DataFrame
    df = pd.DataFrame({"Repeat Count": repeat_counts})
    min_bin, max_bin = min(repeat_counts), max(repeat_counts)

    # Create histogram
    fig = px.histogram(
        df, x="Repeat Count",
        nbins=(max_bin - min_bin + 1),
        title="Allele Size Distribution"
    )

    # Style the layout
    fig.update_layout(
        width=900, height=500,
        margin=dict(t=160),
        title_font_size=18,
        title_font=dict(weight="bold"),
        title_font_color="black",
        xaxis_title="Repeat Count",
        yaxis_title="Allele Count",
        bargap=0.1,
        plot_bgcolor='white'
    )

    # Add x and y axes with zero lines
    fig.update_xaxes(
        zeroline=True,
        zerolinecolor='black'
    )
    fig.update_yaxes(
        zeroline=True,
        zerolinecolor='black'
    )
    
    # Add hover template
    fig.update_traces(
        hovertemplate="Repeat Count=%{x}<br>Allele Count=%{y}<extra></extra>"
    )

    # Get normal and pathogenic ranges from JSON
    benign_min = matched_locus.get("benign_min")
    benign_max = matched_locus.get("benign_max")
    pathogenic_min = matched_locus.get("pathogenic_min")
    pathogenic_max = matched_locus.get("pathogenic_max")

    # Ensure x-axis includes full range
    all_x = repeat_counts + [
        benign_min or 0, benign_max or 0,
        pathogenic_min or 0, pathogenic_max or 0
    ]
    fig.update_layout(xaxis_range=[0, max(all_x) + 10])

    # Calculate data range for pixel-to-unit conversion
    x_range = fig.layout.xaxis.range
    x_span = x_range[1] - x_range[0]

    # Position labels using overlap-aware logic
    label_positions = {}
    if benign_min is not None and benign_max is not None:
        add_range_marker_or_line(
            fig, benign_min, benign_max,
            label_id="normal",
            label_positions=label_positions,
            chart_width_px=fig.layout.width,
            x_span=x_span
        )
    if pathogenic_min is not None and pathogenic_max is not None:
        add_range_marker_or_line(
            fig, pathogenic_min, pathogenic_max,
            label_id="pathogenic",
            label_positions=label_positions,
            chart_width_px=fig.layout.width,
            x_span=x_span
        )

    # Save the plot
    plot_filename = f"{chrom}_{pos}_allele_dist.html"
    fig.write_html(os.path.join(OUTPUT_DIR, plot_filename))
    print(f"Saved plot for {chrom}:{pos}")
