#!/bin/bash/python

file_to_split = "LC1477_18.gb"


def read_file(genbank):
    list1 = []
    list2 = []
    list3 = []
    list4 = []
    round = 0

    with open(genbank, 'r') as file:
        for line in file:
            if not line.startswith("//"):
                if round == 0:
                    list1.append(line)
                elif round == 1:
                    list2.append(line)
                elif round == 2:
                    list3.append(line)
                elif round == 3:
                    list4.append(line)
            else:
                if round == 0:
                    list1.append(line)
                elif round == 1:
                    list2.append(line)
                elif round == 2:
                    list3.append(line)
                elif round == 3:
                    list4.append(line)
                round = round + 1
    dic = {'l1': list1, 'l2': list2, 'l3': list3, 'l4': list4}
    return dic

def write_file(liste, nb):
    if nb == 1:
        genbank = "CP035008.gb"
    elif nb == 2:
        genbank = "CP035009.gb"
    elif nb == 3:
        genbank = "CP035010.gb"
    elif nb == 4:
        genbank = "CP035011.gb"

    with open(genbank, 'w') as file:
        for line in liste:
            file.write(line)

def loop_through_lists(file_to_split):
    lists = read_file(file_to_split)
    nb = 0
    for liste in lists.keys():
        nb = nb + 1
        write_file(lists[liste], nb)

def main(file_to_split):
    loop_through_lists(file_to_split)


main(file_to_split)
