import pandas as pd
import numpy as np
from sqlalchemy import text, create_engine

CLIENT = "postgresql+psycopg2://akshay:8055@127.0.0.1:5432/mitocancer"
engine = create_engine(CLIENT)

data = pd.read_excel(
    "/home/akshay/Projects/MitoCancer/data/Human.MitoCarta3.0.xls",
    sheet_name=1
)

data = data[[
    "HumanGeneID",
    "Symbol",
    "MouseOrthologGeneID",
    "Description",
    "UniProt",
    "MitoCarta3.0_MitoPathways",
    "MitoCarta2.0_Score",
    "TargetP_Score",
    "hg19_Chromosome",
    "EnsemblGeneID_mapping_version_20200130"
]]

data_dict = data.to_dict(orient='records')

def clean(value):
    if value is None:
        return None
    if isinstance(value, float) and np.isnan(value):
        return None
    return value


gene_query = text("""
    INSERT INTO mitocartamaster (
        entrezid,
        genesymbol,
        mouseortholog,
        genedescription,
        uniprotid,
        mitopathway,
        mitocartascore,
        targetpscore,
        hg19chr
    )
    VALUES (
        :entrezid,
        :genesymbol,
        :mouseortholog,
        :description,
        :uniprotid,
        :mitopathway,
        :mitocartascore,
        :targetpscore,
        :hg19chr
    )
    RETURNING geneid
""")

ensembl_query = text("""
    INSERT INTO ensemblidmapping (
        ensemblgeneid,
        geneid,
        idrank
    )
    VALUES (:ensembl_id, :gene_id, :rank)
    ON CONFLICT (ensemblgeneid) DO NOTHING
""")

conflicts = []

with engine.connect() as conn:
    try:
        for elem in data_dict:

            result = conn.execute(gene_query, {
                "entrezid" : clean(elem.get("HumanGeneID")),
                "genesymbol" : clean(elem.get("Symbol")),
                "mouseortholog" : clean(elem.get("MouseOrthologGeneID")),
                "description" : clean(elem.get("Description")),
                "uniprotid" : clean(elem.get("UniProt")),
                "mitopathway" : clean(elem.get("MitoCarta3.0_MitoPathways")),
                "mitocartascore" : clean(elem.get("MitoCarta2.0_Score")),
                "targetpscore" : clean(elem.get("TargetP_Score")),
                "hg19chr" : clean(elem.get("hg19_Chromosome"))
            })

            gene_id = result.fetchone()[0]

            raw_ensembl = elem.get("EnsemblGeneID_mapping_version_20200130")

            if clean(raw_ensembl) is not None:
                ensembl_ids = str(raw_ensembl).split("|")

                for rank, ensembl_id in enumerate(ensembl_ids, start=1):
                    ensembl_id = ensembl_id.strip()
                    if ensembl_id:
                        ensembl_result = conn.execute(ensembl_query, {
                            "ensembl_id": ensembl_id,
                            "gene_id"   : gene_id,
                            "rank"      : rank
                        })

                        if ensembl_result.rowcount == 0:
                            conflicts.append({
                                "ensembl_id": ensembl_id,
                                "gene_id"   : gene_id,
                                "gene_symbol": elem.get("Symbol")
                            })

        conn.commit()
        print(f"Done — {len(data_dict)} genes loaded")
        print(f"{len(conflicts)} Ensembl ID conflicts skipped: {conflicts}")

    except Exception as e:
        conn.rollback()
        print(f"Failed — rolled back. Error: {e}")

