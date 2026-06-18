import pandas as pd
import requests
from sqlalchemy import create_engine, text
import time

CLIENT = "postgresql+psycopg2://akshay:8055@127.0.0.1:5432/mitocancer"
engine = create_engine(CLIENT)
query = "SELECT * FROM mitocartamaster"



def fetch_prot_details(id_dict: dict):

    protein_entry = id_dict.get("entryType", None)
    protein_name = (id_dict.get("proteinDescription", {})
                    .get("recommendedName", {})
                    .get("fullName", {})
                    .get("value", None))
    
    protein_existence = id_dict.get("proteinExistence", None)
    sequence_length = id_dict.get("sequence", {}).get("length", None)
    mol_weight = id_dict.get("sequence", {}).get("molWeight", None)


    return{
        "proteinentry" : protein_entry,
        "proteinname" : protein_name,
        "proteinexistence" : protein_existence,
        "sequncelength" : sequence_length,
        "molweight" : mol_weight
    }


def fetch_is_canonical(id_dict: dict, uniprot_id: str):

    comments_list = id_dict.get("comments", [])

    if not comments_list:
        print(f"No comments sec found for {uniprot_id}")
        return None
    
    alternate_prod = None

    for alt_info in comments_list:
        if alt_info.get("commentType") == "ALTERNATIVE PRODUCTS":
            alternate_prod = alt_info
            return True

    if not alternate_prod:
        print(f"No canonicals for {uniprot_id}")
        return False


def fetch_subcellular_loc(id_dict: dict, uniprot_id: str):

    comments_list = id_dict.get("comments", [])

    if not comments_list:
        print(f"No comments sec found for {uniprot_id}")
        return None

    sub_cellular_block = None

    for sub_info in comments_list:
        if sub_info.get("commentType") == "SUBCELLULAR LOCATION":
            sub_cellular_block = sub_info
            break

    if not sub_cellular_block:
        print(f"No sub cellular info found for {uniprot_id}")
        return None

    sub_cellular_locations = sub_cellular_block.get("subcellularLocations", [])

    try:
        first_location = sub_cellular_locations[0]
        location = first_location.get("location", {}).get("value", None)
        return{
            "uniprousubcellloc" : location
        }
    except(IndexError, KeyError, AttributeError):
        None


def add_disease_details_to_db(id_dict: dict, uniprot_id : str, gene_id : str):

    protein_comments = id_dict.get("comments", [])
    if not protein_comments:
        print(f"No comments for {uniprot_id}")
        return None
    
    disease_dict = []

    for dis in protein_comments:
        if dis.get("commentType", {}) == "DISEASE":
            disease_dict.append(dis)
    if not disease_dict:
        print(f"No disease data for {uniprot_id}")

    results = []

    for disease in disease_dict:
        comment_disease = disease.get("disease", {})
        if not comment_disease:
            continue

        description = comment_disease.get("description", None) if comment_disease else None

        if description is None:
            iscancer = None
        elif any(k in description.lower() for k in ("cancer", "tumour", "tumor")):
            iscancer = True
        else:
            iscancer = False

        dis_query = """
            INSERT INTO uniprotMitoDiseaseData(
                geneid,
                uniprotid,
                diseaseid,
                diseaseacronym,
                mimid,
                diseasedescription,
                iscancer
            ) VALUES (
                :geneid,
                :uniprotid,
                :diseaseid,
                :diseaseacronym,
                :mimid,
                :diseasedescription,
                :iscancer
            )
        """
        try:
            with engine.begin() as conn:
                conn.execute(text(dis_query), {
                "geneid" : gene_id,
                "uniprotid" : uniprot_id,
                "diseaseid" : comment_disease.get("diseaseId", None),
                "diseaseacronym" : comment_disease.get("diseaseAccession", None),
                "mimid" : (comment_disease.get("diseaseCrossReference", {}).get("id", None)
                            if comment_disease.get("diseaseCrossReference", {}).get("database") == "MIM" else None),
                "diseasedescription" : comment_disease.get("description", None),
                "iscancer" : iscancer
                }) 
            results.append(f"Done: {comment_disease.get('diseaseId')}")
        
        except Exception as e:
            print(str(e))
            results.append(f"Failed: {comment_disease.get('diseaseId')}")
        
    return f"{uniprot_id}: " + "; ".join(results)



def main():

    df = pd.read_sql(query, engine)

    for gene_id, uniprot_id in zip(df["geneid"], df["uniprotid"]):
        time.sleep(1)
        print(f"Fetching: {uniprot_id}")
        response = requests.get(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json")
        if response.status_code == 200:

            json_data = response.json()
            protein_info = fetch_prot_details(id_dict=json_data)
            location_info = fetch_subcellular_loc(id_dict=json_data,uniprot_id=uniprot_id)
            

            insert_query = """
                INSERT INTO uniprotMitoData(
                    uniprotid,
                    geneid,
                    proteinentry,
                    proteinname,
                    proteinexistence,
                    subcellularlocation,
                    sequencelength,
                    molweight,
                    iscanonical
                ) VALUES (
                    :uniprotid,
                    :geneid,
                    :proteinentry,
                    :proteinname,
                    :proteinexistence,
                    :subcellularlocation,
                    :sequencelength,
                    :molweight,
                    :iscanonical
                )
            """

            with engine.begin() as conn:
                conn.execute(text(insert_query), {
                    "uniprotid" : uniprot_id,
                    "geneid" : gene_id,
                    "proteinentry" : protein_info.get("proteinentry"),
                    "proteinname": protein_info.get("proteinname"),
                    "proteinexistence" : protein_info.get("proteinexistence"),
                    "subcellularlocation" : location_info.get("uniprousubcellloc") if location_info else None,
                    "sequencelength" : protein_info.get("sequncelength"),
                    "molweight" : protein_info.get("molweight"),
                    "iscanonical" : fetch_is_canonical(id_dict=json_data,uniprot_id=uniprot_id)
                })

            disease_db_feed = add_disease_details_to_db(
                id_dict=json_data,
                uniprot_id=uniprot_id, 
                gene_id= gene_id
            )
            print(disease_db_feed)


        else:
            print(f"Failed for {uniprot_id}: "f"{response.status_code}")
            location_info = None
            protein_info = None


if __name__ == "__main__":
    main()
