import argparse
from pathlib import Path

import decoupler as dc
import numpy as np
import scanpy as sc
from bokeh.palettes import Inferno256
from bokeh.plotting import show

from adata_processing_utils import filter_adata, load_brain_blood_data
from plot_utils import interactive_embedding

msigdb = dc.get_resource("MSigDB")

msigdb_wiki = msigdb[msigdb["collection"] == "wikipathways"]
msigdb_reac = msigdb[msigdb["collection"] == "reactome_pathways"]
msigdb_wiki = msigdb_wiki[~msigdb_wiki.duplicated(["geneset", "genesymbol"])]
msigdb_reac = msigdb_reac[~msigdb_reac.duplicated(["geneset", "genesymbol"])]


def main(project_path, cohort, threshold_targets):

    project_path = Path(project_path)

    blood_adata, brain_adata, blood_adata_hvg, brain_adata_hvg = load_brain_blood_data(
        cohort, project_path
    )

    brain_adata_filt = filter_adata(brain_adata, brain_adata_hvg)
    blood_adata_filt = filter_adata(blood_adata, blood_adata_hvg)

    brain_adata_hvg_filt = filter_adata(brain_adata_hvg, brain_adata_hvg)
    blood_adata_hvg_filt = filter_adata(blood_adata_hvg, blood_adata_hvg)

    for i, msigdbs in zip(["wikipathways", "reactome"], [msigdb_wiki, msigdb_reac]):
        dc.run_ora(
            mat=brain_adata_filt,
            net=msigdbs,
            source="geneset",
            target="genesymbol",
            verbose=True,
            n_up=np.ceil(threshold_targets * brain_adata_filt.shape[1]),
            use_raw=False,
        )
        dc.run_ora(
            mat=blood_adata_filt,
            net=msigdbs,
            source="geneset",
            target="genesymbol",
            n_up=np.ceil(threshold_targets * blood_adata_filt.shape[1]),
            verbose=True,
            use_raw=False,
        )
        brain_adata_hvg_filt.obsm[f"ora_estimate_{i}"] = brain_adata_filt.obsm[
            "ora_estimate"
        ].copy()
        blood_adata_hvg_filt.obsm[f"ora_estimate_{i}"] = blood_adata_filt.obsm[
            "ora_estimate"
        ].copy()

    acts_brs = []
    acts_bls = []
    for msigdbs in ["wikipathways", "reactome"]:

        acts_br = dc.get_acts(
            brain_adata_hvg_filt, obsm_key=f"ora_estimate_{msigdbs}"
        ).copy()
        acts_br_v = acts_br.X.ravel()
        max_br_e = np.nanmax(acts_br_v[np.isfinite(acts_br_v)])
        acts_br.X[~np.isfinite(acts_br.X)] = max_br_e
        acts_brs.append(acts_br.copy())

        acts_bl = dc.get_acts(
            blood_adata_hvg_filt, obsm_key=f"ora_estimate_{msigdbs}"
        ).copy()
        acts_bl_v = acts_bl.X.ravel()
        max_bl_e = np.nanmax(acts_bl_v[np.isfinite(acts_bl_v)])
        acts_bl.X[~np.isfinite(acts_bl.X)] = max_bl_e
        acts_bls.append(acts_bl.copy())

    out_path = project_path / cohort

    wiki_paths = [
            "WP_CELLS_AND_MOLECULES_INVOLVED_IN_LOCAL_ACUTE_INFLAMMATORY_RESPONSE",
            "WP_COMPLEMENT_SYSTEM_IN_NEURONAL_DEVELOPMENT_AND_PLASTICITY",
            "WP_TYROBP_CAUSAL_NETWORK_IN_MICROGLIA",
            "cell_type",
        ]

    reac_paths = ["REACTOME_CD163_MEDIATING_AN_ANTI_INFLAMMATORY_RESPONSE", 
                "REACTOME_ANTI_INFLAMMATORY_RESPONSE_FAVOURING_LEISHMANIA_PARASITE_INFECTION",
            "REACTOME_LEISHMANIA_INFECTION",
            "cell_type"]

    for wiki_path in wiki_paths[:-1]:
        acts_brs[0].obs[wiki_path]=acts_brs[0].obsm['ora_estimate_wikipathways'][wiki_path]
        acts_bls[0].obs[wiki_path]=acts_bls[0].obsm['ora_estimate_wikipathways'][wiki_path]
    
    for reac_path in reac_paths[:-1]:
        acts_brs[1].obs[reac_path]=acts_brs[1].obsm['ora_estimate_reactome'][reac_path]
        acts_bls[1].obs[reac_path]=acts_bls[1].obsm['ora_estimate_reactome'][reac_path]

    for wiki_path in wiki_paths[:-1]:
        p = interactive_embedding(acts_brs[0], label=wiki_path,palette_cont=Inferno256, title_font_size="20px")
        show(p)
        p2 = interactive_embedding(acts_bls[0], label=wiki_path,palette_cont=Inferno256, title_font_size="20px")
        show(p2)
    for reac_path in reac_paths[:-1]:
        p = interactive_embedding(acts_brs[1], label=reac_path, palette_cont=Inferno256, title_font_size="20px")
        show(p)
        p2 = interactive_embedding(acts_bls[1], label=reac_path, palette_cont=Inferno256, title_font_size="20px")
        show(p2)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p",
        "--project_path",
        help="Path were data are stored",
        default="/sc/arion/projects/psychgen/lbp/data/ScProcesses_brainBlood_nicole/",
    )
    parser.add_argument(
        "-b", "--cohort", help="cohort identifier", choices=["cohort_1", "cohort_2"]
    )

    parser.add_argument(
        "-t", "--threshold_targets", help="decoupler ORA threshold", type=float, default=0.05
    )

    args = parser.parse_args()
    main(
        args.project_path,
        args.cohort,
        args.threshold_targets,
    )
