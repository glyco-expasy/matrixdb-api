import re


class NetworkManager:

    def __init__(self, database_manager, cache_manager, protein_data_manager, interaction_data_manager=None):
        # initialize database connection
        self.database_manager = database_manager
        self.meta_data_cache = cache_manager.get_meta_data()
        self.protein_data_manager = protein_data_manager
        self.interaction_data_manager = interaction_data_manager

    def merge_associations(self, *interactions):
        result = interactions[0].copy()

        for d in interactions[1:]:
            for key, value in d.items():
                if key in result and isinstance(result[key], list) and isinstance(value, list):
                    result[key] += value
                elif key in result and isinstance(result[key], dict) and isinstance(value, dict):
                    result[key] = self.merge_associations(result[key], value)
                else:
                    result[key] = value

        return result

    def get_interactors(self, interactor_ids, interactor_mapping):
        core_database_connection = self.database_manager.get_primary_connection()
        # Extract the interactors
        print("reading interactors")
        interactors = list(core_database_connection["biomolecules"].find({
            "id": {"$in": list(interactor_ids)}
        }))

        # Check for isoforms (*-1) and pro peptides (*-pro)
        proteoform_interactors = list()
        for interactor_id in interactor_ids:
            if "-PRO" in interactor_id:
                protein_id = interactor_id.split("-PRO")[0]
                protein = core_database_connection["biomolecules"].find_one({
                    'id': protein_id
                })
                protein["id"] = interactor_id
                proteoform_interactors.append(protein)
                #proteoform_interactors.append({
                #    "id": interactor_id,
                #    "name": protein[""]
                #    "type": "protein",
                #    "proteoform": "pro-peptide"
                #})

            elif "-" in interactor_id:
                protein_id = interactor_id.split("-")[0]
                protein = core_database_connection["biomolecules"].find_one({
                    'id': protein_id
                })
                protein["id"] = interactor_id
                proteoform_interactors.append(protein)
                #proteoform_interactors.append({
                #    "id": interactor_id,
                #    "type": "protein",
                #    "proteoform": "isoform"
                #})

        interactors.extend(proteoform_interactors)

        print("Interactors read")
        transformed_interactors = self.transform_interactors(interactors, interactor_mapping)

        # Normalize the interactors
        print("interactors are done")
        return transformed_interactors

    def get_interactions(self, interaction_ids, interactor_mapping):
        core_database_connection = self.database_manager.get_primary_connection()
        # Extract the interactions

        # Extracts a map of experiments
        experiments = dict()
        experiment_map = dict()
        type_map = {
            "1": "Only Experimental",
            "2": "Only Predicted",
            "3": "Experimental and Predicted",
        }
        transformed_interactions = list()
        psimi_names = dict()
        psimi_name_mapping = dict()
        for interaction in core_database_connection["interactions"].aggregate([
                {
                    "$match": {
                      "id": {"$in": list(interaction_ids)}
                    }
                },
                {
                    "$lookup": {
                        "from": "experiments",
                        "localField": "experiments.direct.binary",
                        "foreignField": "id",
                        "pipeline": [{"$project": {"interaction_detection_method": 1}}],
                        "as": "binary_detection_method"
                    }
                },
                {
                    "$lookup": {
                        "from": "experiments",
                        "localField": "experiments.direct.spoke_expanded_from",
                        "foreignField": "id",
                        "pipeline": [{"$project": {"interaction_detection_method": 1}}],
                        "as": "spoke_detection_method"
                    }
                }
        ]):
            transformed_interaction = {
                "id": interaction["id"],
                "participants": [interactor_mapping[participant] for participant in interaction["participants"] if participant in interactor_mapping],
                "experiments": {
                    "direct": {
                        "binary": list(),
                        "spoke_expanded_from": list()
                    }
                }
            }
            if "experiments" in interaction:
                if "binary" in interaction["experiments"]["direct"]:
                    for binary_experiment in interaction["experiments"]["direct"]["binary"]:
                        if binary_experiment not in experiments:
                            experiments[binary_experiment] = len(experiments) + 1
                            experiment_map[experiments[binary_experiment]] = binary_experiment
                        transformed_interaction["experiments"]["direct"]["binary"].append(experiments[binary_experiment])

                if "spoke_expanded_from" in interaction["experiments"]["direct"]:
                    for spoke_exapnded in interaction["experiments"]["direct"]["spoke_expanded_from"]:
                        if spoke_exapnded not in experiments:
                            experiments[spoke_exapnded] = len(experiments) + 1
                            experiment_map[experiments[spoke_exapnded]] = spoke_exapnded
                        transformed_interaction["experiments"]["direct"]["spoke_expanded_from"].append(experiments[spoke_exapnded])

            if 'score' in interaction:
                if interaction['score'] != '-':
                    transformed_interaction['score'] = interaction['score']

            if 'prediction' in interaction and 'experiments' in interaction:
                transformed_interaction['type'] = 3

            if 'prediction' in interaction and 'experiments' not in interaction:
                transformed_interaction['type'] = 2

            if 'prediction' not in interaction and 'experiments' in interaction:
                transformed_interaction['type'] = 1

            interaction_detection_methods = list()
            if "binary_detection_method" in interaction:
                if len(interaction["binary_detection_method"]) > 0:
                    for detection_method in interaction["binary_detection_method"]:
                        # Extract the psimi id, should be fixed in data
                        if "MI" in detection_method["interaction_detection_method"]:
                            match = re.search(r'MI:\d+', detection_method["interaction_detection_method"])
                            if match:
                                mi_id = match.group()
                        else:
                            mi_id = detection_method["interaction_detection_method"]

                        psimi_name = self.meta_data_cache["psimi"][mi_id]["name"]
                        if psimi_name not in psimi_names:
                            psimi_names[psimi_name] = len(psimi_names) + 1
                            psimi_name_mapping[psimi_names[psimi_name]] = psimi_name

                        interaction_detection_methods.append(psimi_names[psimi_name])

            if "spoke_detection_method" in interaction:
                if len(interaction["spoke_detection_method"]) > 0:
                    for detection_method in interaction["spoke_detection_method"]:
                        # Extract the psimi id, should be fixed in data
                        if "psi-mi" in detection_method["interaction_detection_method"]:
                            match = re.search(r'MI:\d+', detection_method["interaction_detection_method"])
                            if match:
                                mi_id = match.group()
                        else:
                            mi_id = detection_method["interaction_detection_method"]

                        psimi_name = self.meta_data_cache["psimi"][mi_id]["name"]
                        if psimi_name not in psimi_names:
                            psimi_names[psimi_name] = len(psimi_names) + 1
                            psimi_name_mapping[psimi_names[psimi_name]] = psimi_name

                        interaction_detection_methods.append(psimi_names[psimi_name])
            transformed_interaction["detection_method"] = interaction_detection_methods

            transformed_interactions.append(transformed_interaction)

        return {
            "interactions": transformed_interactions,
            "context": {
                "interactions": {
                    "type": type_map,
                    "detection_method": psimi_name_mapping,
                },
                "experiment_mapping": experiment_map
            }
        }

    def transform_interactors(self, interactors, interactor_mapping):
        # Extract the expressions for interactors
        gene_expressions = self.protein_data_manager.get_gene_expression_for_proteins(
            list(i['id'] for i in interactors))
        print("gene expression read")
        proteomics_expressions = self.protein_data_manager.get_proteomics_expressions_for_proteins(
            list(i['id'] for i in interactors))
        print("proteomics expression are read")
        transformed_interactors = list()
        tissue_names = dict()
        tissue_name_mapping = dict()

        proteomics_expression_tissue_names = dict()
        prot_expression_tissue_name_mapping = dict()

        proteomics_expression_sample_names = dict()
        prot_expression_sample_name_mapping = dict()

        for interactor in interactors:

            transformed_interactor = {
                'id': interactor['id'],
                'type': interactor['type']
            }

            if 'proteoform' in interactor:
                transformed_interactor['proteoform'] = interactor['proteoform']

            if 'names' in interactor:
                if 'recommended_name' in interactor['names']:
                    transformed_interactor['name'] = interactor['names']['recommended_name']
                elif 'name' in interactor['names']:
                    transformed_interactor['name'] = interactor['names']['name']

            if 'molecular_details' in interactor and 'pdb' in interactor['molecular_details']:
                transformed_interactor['pdb'] = len(interactor['molecular_details']) > 0

            if 'ecm' in interactor:
                transformed_interactor['ecm'] = 'Yes'

            #if 'relations' in interactor and 'gene_name' in interactor['relations']:
            #    if interactor['relations']['gene_name'] not in gene_names:
            #        gene_names[interactor['relations']['gene_name']] = len(gene_names) + 1
            #    transformed_interactor['geneName'] = gene_names[interactor['relations']['gene_name']]

            # Get expressions and add to interactor
            if interactor['type'] == 'protein':
                if interactor['id'] in gene_expressions:
                    expressions = gene_expressions[interactor['id']]
                    transformed_interactor['geneExpression'] = list()
                    for expression in expressions:
                        tissue_name = expression["tissue"]
                        if tissue_name not in tissue_names:
                            tissue_names[tissue_name] = len(tissue_names) + 1
                            tissue_name_mapping[tissue_names[tissue_name]] = tissue_name

                        transformed_interactor['geneExpression'].append({
                            'tissue': tissue_names[tissue_name],
                            'tpm': expression['tpm']
                        })

                transformed_interactor['proteomicsExpression'] = list()
                if interactor['id'] in proteomics_expressions:
                    expressions = proteomics_expressions[interactor['id']]
                    for expression in expressions:
                        if expression["tissue"] not in proteomics_expression_tissue_names:
                            number_mapping = len(proteomics_expression_tissue_names) + 1
                            proteomics_expression_tissue_names[expression["tissue"]] = number_mapping
                            prot_expression_tissue_name_mapping[number_mapping] = expression["tissue"]

                        if expression["sampleName"] not in proteomics_expression_sample_names:
                            number_mapping = len(proteomics_expression_sample_names) + 1
                            proteomics_expression_sample_names[expression["sampleName"]] = number_mapping
                            prot_expression_sample_name_mapping[number_mapping] = expression["sampleName"]

                        transformed_interactor['proteomicsExpression'].append({
                            'tissue': proteomics_expression_tissue_names[expression["tissue"]],
                            'sampleName': proteomics_expression_sample_names[expression["sampleName"]],
                            'score': expression["score"],
                        })

            transformed_interactor["id"] = interactor_mapping[interactor['id']]
            transformed_interactors.append(transformed_interactor)
        print("done transforming")

        context = {
            "interactors": {
                "interactor_mapping": {value: key for key, value in interactor_mapping.items()},
                "geneExpression": {
                    "tissue": tissue_name_mapping,
                },
                "proteomicsExpression": {
                    "tissue": prot_expression_tissue_name_mapping,
                    "sampleName": prot_expression_sample_name_mapping,
                }
            }
        }

        return {
            "interactors": transformed_interactors,
            "context": context
        }

    def generate_network(self, biomolecules, add_second_neighborhood=True):
        '''
        Generates the interaction network of one or more biomolecules
        Network consists of all biomolecule partners of the give biomolecule list, found in interactions
        :param biomolecules:
        :param add_second_neighborhood:
        :return: {
                    biomolecules: list of all biomoleules
                    interactions: list of all interactions
                 }
        '''
        interactor_mapping = dict()

        # Get the 1-neighborhood of each biomolecule and derive the interactions and partners
        interactor_ids = set()
        interaction_ids = set()

        biomolecules_list = list()
        partner_list = None
        for biomolecule in biomolecules:
            # Add current biomolecule
            interactor_ids.add(biomolecule)

            if biomolecule not in interactor_mapping:
                interactor_mapping[biomolecule] = len(interactor_mapping) + 1

            biomolecules_list.append(interactor_mapping[biomolecule])

            # Partners
            neighborhood = self.interaction_data_manager.get_neighborhood(biomolecule)
            for biomolecule_form in neighborhood:
                partner_list = neighborhood[biomolecule_form]
                interactor_ids.update(partner_list)
                interactor_ids.add(biomolecule_form)

            # First neighborhood interactions
            for biomolecule_form in neighborhood:
                partner_list = neighborhood[biomolecule_form]
                for partner in partner_list:
                    if partner not in interactor_mapping:
                        interactor_mapping[partner] = len(interactor_mapping) + 1

                    sorted_partners = sorted([biomolecule_form, partner])
                    interaction_id = f"{sorted_partners[0]}__{sorted_partners[1]}"
                    interaction_ids.add(interaction_id)

                if biomolecule_form not in interactor_mapping:
                    interactor_mapping[biomolecule_form] = len(interactor_mapping) + 1

                biomolecules_list.append(interactor_mapping[biomolecule_form])

            # Check for self interactions
            if partner_list != None and biomolecule in partner_list:
                interaction_ids.add(f"{biomolecule}__{biomolecule}")

            print(f"First neighborhood: {len(interaction_ids)}")

            # If second neighrohood requested
            second_neighborhood = set()
            if add_second_neighborhood:
                second_neigborhood_interactions = set()
                for interactor_id in interactor_ids:
                    # Build the network adjacency
                    if interactor_id not in interactor_mapping:
                        interactor_mapping[interactor_id] = len(interactor_mapping) + 1

                    neighborhood = self.interaction_data_manager.get_neighborhood(interactor_id)
                    for biomolecule_form in neighborhood:
                        for second_neighbor in neighborhood[biomolecule_form]:
                            # Only consider the interactions among first neighrhoods
                            if second_neighbor not in interactor_ids:
                                continue
                            second_neighborhood.add(second_neighbor)

                            # Second neighborhood interactions
                            sorted_partners = sorted([interactor_id, second_neighbor])
                            interaction_id = f"{sorted_partners[0]}__{sorted_partners[1]}"
                            second_neigborhood_interactions.add(interaction_id)

                            if second_neighbor not in interactor_mapping:
                                interactor_mapping[second_neighbor] = len(interactor_mapping) + 1

                print(f"Second neigbhrhood: {len(second_neighborhood)}")
                print(f"Second neighborhood interaction {len(second_neigborhood_interactions)}")
                interactor_ids.update(second_neighborhood)
                interaction_ids.update(second_neigborhood_interactions)

        interactor_list = self.get_interactors(interactor_ids, interactor_mapping)
        interaction_list = self.get_interactions(interaction_ids, interactor_mapping)

        context = {
            "interactors": {value: key for key, value in interactor_mapping.items()}
        }
        context.update(interactor_list["context"])
        context.update(interaction_list["context"])
        print("Done generating network")

        # Participant count is calculated with following rule
        # take all completely external partners, note that, there could be multiple interactions with same parnter
        # due to isoforms in proteins  a-b, a-b1, but we only count it once
        # for all self interactions with or without isoforms will also count as one parnter
        partner_count = set()
        biomolecule_ids = list(context["interactors"]["interactor_mapping"][biomol] for biomol in biomolecules_list)
        for interaction in interaction_list["interactions"]:
            p1 = interaction['id'].split('__')[0]
            p2 = interaction['id'].split('__')[1]

            if p1 not in biomolecule_ids:
                partner_count.add(p1)

            if p2 not in biomolecule_ids:
                partner_count.add(p2)

            if p1 in biomolecule_ids or p2 in biomolecule_ids:
                if '-' in p1:
                    p1 = p1.split('-')[0]

                if '-' in p2:
                    p2 = p2.split('-')[0]

                if p1 == p2:
                    partner_count.add(p1)

        return {
            "partnerCount": len(partner_count),
            "interactors": interactor_list["interactors"],
            "interactions": interaction_list["interactions"],
            "biomolecules": biomolecules_list,
            "context": context,
        }


