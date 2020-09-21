from genericpath import isfile
from os import listdir
from os.path import join

from Bio.KEGG.REST import kegg_get

#from hw1 import OUTPUT_DIR


def get_KEGG_data(database, query, option=None):

    # first check to see if we have downloaded the query into a file
    all_files_in_output = [f for f in listdir("output") if isfile(join("output", f))]
    if "{}.txt".format(query) not in all_files_in_output:
        http_response = kegg_get(query, option)
        response_string  = http_response.read()

        with open("{}/{}.txt".format("output", query), "w") as file_to_save:
            file_to_save.write(response_string)

    file_relative_path_and_name = "{}/{}.txt".format("output", query)
    # with open(file_relative_path_and_name, "r") as queried_file:
    #     for line in queried_file.readlines():
    #         print line

    return file_relative_path_and_name

