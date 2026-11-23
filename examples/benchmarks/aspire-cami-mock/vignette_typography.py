"""Presentation-only typography for the standalone CAMI vignette SVGs."""
import re


def style_svg(svg, name):
    # Matplotlib groups each title (including multiline titles) separately.
    # Match its displayed text so that ticks, legends and numeric annotations survive.
    titles = {
        "Bray_Boxplot_sample_type_within_patient": ["Within-Patient Bray-Curtis by Sample-Type Pair"],
        "Distance_contrasts_ABC": ["Bray-Curtis Distances Across Planned Contrasts (A/B/C)"],
        "alpha_faceted_Shannon_Case": ["Type_Group = Airways", "Type_Group = Oral"],
        "community": ["Within-library diversity", "Bray–Curtis ordination"],
        "final_micro_upset": ["ASV Membership by Type_Group (Final)"],
        "power_curves_cancer_vs_control_simple": ["Power Analysis: Cancer vs Control (Observed Effect)", "Airways", "Oral"],
    }

    def clean_group(match):
        group = match.group()
        labels = re.findall(r"<text\b[^>]*>([^<]+)</text>", group)
        if not labels:
            return group
        label = labels[0]
        is_title = (
            label in titles.get(name, ())
            or (name.startswith("network") and label.startswith("SPIEC-EASI Network"))
            or (name == "site-case-networks" and " network indicators" in label)
            or (name == "case-isa-single-group" and label.startswith("Case vs control ·"))
            or (name == "patient-isa-acetone" and label.startswith("Airway indicator and acetone across "))
        )
        # The caption identifies the same alpha reference line in all four panels.
        if name == "power_curves_cancer_vs_control_simple" and label == "α=0.05":
            return ""
        return "" if is_title else group

    svg = re.sub(r'(?ms)^([ ]*)<g id="text_\d+"[^>]*>.*?^\1</g>\n?', clean_group, svg)
    if name == "mock-read-flow":
        svg = re.sub(r"<text\b[^>]*>(?:Data Loss Flow|Read pairs / merged-read equivalents; band widths represent counts\.)</text>", "", svg)
    # Large network plates need larger legend type at the same page width.
    network_key = name.startswith(("network_", "network-degree"))

    def font(match):
        original = float(match[1])
        if network_key and original >= 17:
            return "font-size: 17px"
        size = 12 if original >= 12 else 10
        return f"font-size: {size}px"

    return re.sub(r"font-size:\s*([\d.]+)px", font, svg)
