import os

from model import Reaction, RClass
from model.RClass import parse, Record, RClass_cross_reference
from network.rest import get_KEGG_data

ROOT_DIR = ""
OUTPUT_DIR = ""

def setup_env():
    global ROOT_DIR
    global  OUTPUT_DIR
    ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
    OUTPUT_DIR = ROOT_DIR + "/output"

    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

def test_get_data():
    get_KEGG_data("cpd", "C00154")
    get_KEGG_data("ec", "1.2.1.50")
    get_KEGG_data("rc", "RC00184")
    get_KEGG_data("cpd", "C00150")
    get_KEGG_data("rn", "R01277")


def test_part2_compound():
    with open(OUTPUT_DIR + "/C00154.txt") as handle:
         for record in parse(handle):
            print("%s %s" % (record.entry, record.name[0]))

def test_part2_reaction():
    with open(OUTPUT_DIR + "/R01277.txt") as handle:
         for record in parse(handle):
            # print("%s %s %s %s" % (record.entry, record.name[0], record.enzyme, record.definition))
            print record

def test_part3_rclass():
    with open(OUTPUT_DIR + "/RC00099.txt") as handle:
         for record in parse(handle):
            # print("%s %s %s %s" % (record.entry, record.name[0], record.enzyme, record.definition))
            print record.reaction
            print record.rpair
            print record

# this is just to test our part 4 results
def test_part4_rclass_part5_find_common_rpairs():
    record1 = Record()
    record1 = record1.analyze_RClass("RC00184")
    record2 = Record()
    record2 = record2.analyze_RClass("RC00300")
    record3 = Record()
    record3 = record3.analyze_RClass("RC00099")

    rpair = [record1.rpair, record2.rpair, record3.rpair]

    for r in rpair:
        print r
        print len(r)

    rpair = [["1"], ["1"], ["2"]]

    common_0_1 = set()
    common_0_2 = set()
    common_1_2 = set()
    common_0_1_2 = set()

    for s_p in rpair[0]:
        if s_p in rpair[1]:
            common_0_1.add(s_p)
            if s_p in rpair[2]:
                common_0_1_2.add(s_p)
                common_0_2.add(s_p)
        elif s_p in rpair[2]:
            common_0_2.add(s_p)
    for s_p in rpair[1]:
        if s_p in rpair[2]:
            common_1_2.add(s_p)

    print("common 0 and 1: ", common_0_1)
    print("common 0 and 2: ", common_0_2)
    print("common 1 and 2: ", common_1_2)
    print("common 0 and 1 and 2: ", common_0_1_2)


INPUT_RCLASSES = ["RC00184", "RC00300", "RC00099"]

def test_part5():

    for rc_number in INPUT_RCLASSES:
        cross_reference_rc(rc_number)




def cross_reference_rc(rc_number):
    print("RClass {}").format(rc_number)
    file_name = get_KEGG_data("rc", rc_number)
    cr_ds_list = []
    with open(file_name) as handle:
        for rclass_record in RClass.parse(handle):
            for rn in rclass_record.reaction:
                file_name_for_reaction = get_KEGG_data("rn", rn)
                with open(file_name_for_reaction) as handle_reaction:
                    for reaction_record in Reaction.parse(handle_reaction):
                        # cross-reference data structure
                        cr_ds = RClass_cross_reference(reaction_record, rclass_record.entry)
                        cr_ds_list.append(cr_ds)
    output_tab_separated_file_path_and_name = OUTPUT_DIR + '/results_tab_delimited.txt'
    output_csv_file_path_and_name = output_tab_separated_file_path_and_name.replace("_tab_delimited.txt", ".csv")
    with open(output_tab_separated_file_path_and_name, 'w') as fp:

        fp.write("reaction\tenzymes\trpairs\tk-number\tk-names\n")
        for item in cr_ds_list:
            fp.write(str(item))

    # To use pretty table for nice output. Need to convert to CSV first.
    tab_csv(output_csv_file_path_and_name, output_tab_separated_file_path_and_name)
    fix_pretty_table_csv_bug(output_csv_file_path_and_name)
    from prettytable import from_csv
    fp = open(output_csv_file_path_and_name, "r")
    mytable = from_csv(fp)
    fp.close()
    print mytable, "\n"


def fix_pretty_table_csv_bug(output_csv_file_path_and_name):
    """
    BUG 1:
    prettytable library has an odd bug where occurance of "- " in CSVs would crash the program, therefore we manually
    replace the mentioned expression with " -".

    BUG2:
    Last CVS line should not be empty

    :param output_csv_file_path_and_name: path to the csv file to be fixed
    :return: nothing
    """
    csv_lines_output = []
    fp_read = open(output_csv_file_path_and_name, "r")
    csv_lines = fp_read.readlines()
    for csv_line in csv_lines:
        csv_line = csv_line.replace("- ", " -")
        csv_lines_output.append(csv_line)
    fp_read.close()
    fp_write = open(output_csv_file_path_and_name, "w")
    for csv_line in csv_lines_output:
        fp_write.write(csv_line)
    fp_write.close()


def tab_csv(output_csv_file_path_and_name, output_tab_separated_file_path_and_name):
    import csv
    txt_file = output_tab_separated_file_path_and_name
    csv_file = output_csv_file_path_and_name
    in_txt = csv.reader(open(txt_file, "rb"), delimiter="\t")
    out_csv = csv.writer(open(csv_file, "wb"))
    out_csv.writerows(in_txt)


if __name__ == "__main__":

    setup_env()

    test_part5()




