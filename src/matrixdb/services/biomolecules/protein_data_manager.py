from itertools import groupby

from src.matrixdb.utils.database.database_manager import DatabaseManager



def tissue_maxima(expression_values):
    """One entry per tissue: the sample with the highest TPM.

    `tpm` is stored as a string, so it must be compared as a number: as strings "9.5" is greater than "123.4", which
    made the API return the wrong sample for about 40% of the tissue values. Tissues are ordered case-insensitively so the
    order is the same as the one the UI uses for its labels.
    """
    ordered = sorted(expression_values, key=lambda x: x["uberonName"].strip('"').lower())
    result = list()
    for tissue_uberon_name, group in groupby(ordered, key=lambda x: x["uberonName"]):
        max_group = max(group, key=lambda x: float(x["tpm"]))
        result.append({
            "tissueId": max_group["tissueId"],
            "tissue": tissue_uberon_name.strip('"'),
            "tpm": max_group["tpm"]
        })
    return result

class ProteinDataManager:

    def __init__(self, app_config, cache_manager, database_manager):
        self.app_config = app_config
        self.meta_data_cache = cache_manager.get_meta_data()
        self.database_manager = database_manager

    def get_gene_expressions(self, protein_id):
        secondary_database_connection = self.database_manager.get_secondary_connection()
        expressions_by_protein = secondary_database_connection["geneExpression"].find_one({
            "uniprot": protein_id
        })

        gene_expressions = list()
        if expressions_by_protein is not None:

            expression_values = list()
            for expression in expressions_by_protein["expressions"]:
                for tissue_id, expression_value in expression.items():
                    expression_by_tissue = expression_value[0]
                    expression_by_tissue['tissueId'] = tissue_id
                    expression_values.append(expression_by_tissue)

            gene_expressions = tissue_maxima(expression_values)
        return gene_expressions

    def get_gene_expression_for_proteins(self, protein_ids):
        secondary_database_connection = self.database_manager.get_secondary_connection()
        expressions_by_proteins = secondary_database_connection["geneExpression"].find({
            "uniprot": {
                '$in': protein_ids
            }
        })

        gene_expressions = dict()
        if expressions_by_proteins is not None:
            for expressions_by_protein in expressions_by_proteins:
                protein = expressions_by_protein['uniprot']
                gene_expressions[protein] = list()
                expression_values = list()
                for expression in expressions_by_protein["expressions"]:
                    for tissue_id, expression_value in expression.items():
                        expression_by_tissue = expression_value[0]
                        expression_by_tissue['tissueId'] = tissue_id
                        expression_values.append(expression_by_tissue)

                gene_expressions[protein] = tissue_maxima(expression_values)
        return gene_expressions

    def get_proteomics_expressions(self, protein_id):
        secondary_database_connection = self.database_manager.get_secondary_connection()
        proteomics_expression_by_protein = secondary_database_connection["proteomicsExpressions"].find_one({
            "uniprot": protein_id
        })

        prot_expressions = list()
        if proteomics_expression_by_protein is not None:
            for expression in proteomics_expression_by_protein["expressions"]:
                tissue_name = self.get_tissue_name(expression)
                if tissue_name is not None:
                    prot_expressions.append(
                        {
                            "tissue": tissue_name,
                            "sampleName": expression["sampleName"] if "sampleName" in expression else "",
                            "sample": expression["sample"] if "sample" in expression else "",
                            "score": expression["nsaf"] if "nsaf" in expression else 0
                        })

        return prot_expressions

    def get_proteomics_expressions_for_proteins(self, protein_ids):
        secondary_database_connection = self.database_manager.get_secondary_connection()
        proteomics_expression_by_proteins = secondary_database_connection["proteomicsExpressions"].find({
            "uniprot": {
                '$in': protein_ids
            }
        })

        prot_expressions = dict()
        if proteomics_expression_by_proteins is not None:
            for proteomics_expression_by_protein in list(proteomics_expression_by_proteins):
                prot_expressions[proteomics_expression_by_protein['uniprot']] = list()
                for expression in proteomics_expression_by_protein["expressions"]:

                    tissue_name = self.get_tissue_name(expression)
                    if tissue_name is not None:
                        prot_expressions[proteomics_expression_by_protein['uniprot']].append(
                            {
                                "tissue": tissue_name,
                                "sampleName": expression["sampleName"] if "sampleName" in expression else "",
                                "sample": expression["sample"] if "sample" in expression else "",
                                "score": expression["nsaf"] if "nsaf" in expression else 0,
                            })

        return prot_expressions

    def get_tissue_name(self, expression_data):
        if "BTO" in expression_data["tissueId"]:
            if expression_data['tissueId'] in self.meta_data_cache["bto"]:
                return self.meta_data_cache["bto"][expression_data["tissueId"]]["name"]

        if "UBERON" in expression_data["tissueId"]:
            if expression_data['tissueId'] in self.meta_data_cache["uberon"]:
                return self.meta_data_cache["uberon"][expression_data["tissueId"]]["name"]

    def get_gene_name(self, protein_id):
        primary_connection = self.database_manager.get_primary_connection()
        protein = primary_connection['biomolecules'].find_one({
            'id': protein_id
        })
        if type(protein['relations']['gene_name']) is list:
            return protein['relations']['gene_name'][0]
        else:
            return protein['relations']['gene_name']