from simplesam import Reader
import argparse
from collections import defaultdict
import re
import sys

class Read:
  def __init__(self, my_read):
    self.chr=my_read.rname
    self.aligned_ref=''.join(my_read.parse_md())
    self.seq=my_read.gapped('seq')
    aligned_coords=my_read.coords
    self.start_coord=aligned_coords[0]-1 # we want the 0-BASED genomic coordinates
    self.end_coord=self.start_coord+len(self.aligned_ref)-1
    self.methylation_string=my_read.gapped('XM')
    self.read_name=my_read.qname

# code from https://stackoverflow.com/questions/4836710/is-there-a-built-in-function-for-string-natural-sort
def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)

# def find_target_cgs(mychr,mystart,mystop,cg_dict):
class Target:
    order=''

    def __init__(self,chr,pos,probe):
        self.chr=chr
        self.pos=int(pos)
        self.probe=probe
        self.name=chr+':'+str(pos)
    
    def set_order(self,order):
        self.order=order

def get_methylation(my_read,my_targets):
    methylation={}

    # are any targets found in this read?
    for target in my_targets:
        if my_targets[target].chr==my_read.chr and my_targets[target].pos > my_read.start_coord and my_targets[target].pos <= my_read.end_coord:
            meth_index=my_targets[target].pos-my_read.start_coord-1
            methylation_state=my_read.methylation_string[meth_index]

            # print('cg: {}; start: {}, stop: {}, meth: {}, seq: {}, cg_meth: {}, cg_base: {}'.format(my_targets[target].name,my_read.start_coord, my_read.end_coord,my_read.methylation_string,my_read.seq,my_read.methylation_string[meth_index],my_read.seq[meth_index]))
            # print('\n')

            if my_read.aligned_ref[meth_index:(meth_index+2)].upper()=='CG':
                # Z=1; z=0
                if methylation_state=='Z':
                    methylation[my_targets[target].name]='1'
                elif methylation_state=='z':
                    methylation[my_targets[target].name]='0'
            else:
                print("ERROR...target CG position for {} does not match a CG site".format(my_targets[target].name))
                print(methylation_state)
                print(meth_index)
                print(my_read.methylation_string)
                print(my_read.seq)
                print(my_read.aligned_ref)
                print(my_read.chr)
                print(my_read.start_coord)
                print(my_targets[target].pos)
                print(my_read.end_coord)
                print("EXITING...")
                sys.exit()

    return(methylation)






if __name__ == "__main__" :

    ###############
    ## ARGUMENTS ##
    ###############
    parser=argparse.ArgumentParser(description='Find methylation on same paired-end illumina sequencing molecules')
    parser.add_argument('--sam', help='.sam file',type=str, required=True)
    parser.add_argument('--targets', help='target file',type=str, required=True)

    args=parser.parse_args()

    targets_file_name=args.targets
    sam_file_name=args.sam

    output_file_name=sam_file_name+'-methylation_results.tsv'

    ##########################################
    ## READ TARGET CGs / FORMAT OUTPUT FILE ##
    ##########################################
    targets_dict=defaultdict(list)
    target_objects={}
    # target_to_probe_dict={}
    with open(targets_file_name,'r') as targets_file:
        targets_file.readline()

        for line in targets_file:
            line=line.strip('\n').split('\t')
            chr=line[0]
            start=line[1]
            probe=line[3]
            name=chr+':'+start
            
            curr_target=Target(chr,start,probe)
            target_objects[name]=curr_target

            targets_dict[chr].append(start)
            targets_dict[chr]=sorted(targets_dict[chr]) # make sure CGs are sorted; this is so CGs in output will be sorted by coordinates

    # print(target_objects)
    ## FORMAT OUTPUT FILE
    # get list of sorted chr
    sorted_chr=natural_sort(targets_dict.keys())

    output_file=open(output_file_name,'w')
    header_line=['Read_Name']
    header_order={}
    
    counter=1
    for chr in sorted_chr:
        for pos in targets_dict[chr]:
            cg_name=str(chr)+':'+str(pos)
            header_line.append(cg_name)
            
            # add order to target objects
            target_objects[cg_name].set_order(counter)
                
            counter+=1

    # print(header_line)
    # print(header_order)

    # for target in target_objects:
    #     print(target_objects[target].name,target_objects[target].chr,target_objects[target].pos,target_objects[target].probe,target_objects[target].order)

    output_file.write('\t'.join(header_line)+'\n')

    #######################
    ## ANALYZE EACH READ ##
    #######################
    sam_file=Reader(open(sam_file_name, 'r'))

    for forward_read in sam_file:
        if forward_read.mapped and forward_read.paired: # makes sure reads map to a reference seq and are paired
                reverse_read=sam_file.next()

                if forward_read.qname != reverse_read.qname:
                    print("READS NOT PAIRED, out of sync")
                    print("forward_read: ",forward_read.qname,"reverse_read: ",reverse_read.qname)
                    print("EXITING...")
                    sys.exit()


                # print('forward read: {}'.format(forward_read))
                forward_read=Read(forward_read)

                # print('chr: {}, start coord: {}, end coord: {}, meth: {}, aligned seq: {}, refseq: {}'.format(forward_read.chr,forward_read.start_coord,forward_read.end_coord,forward_read.methylation_string,forward_read.aligned_ref,forward_read.seq))
                

                # print('reverse read: {}'.format(reverse_read))
                reverse_read=Read(reverse_read)
                # print('chr: {}, start coord: {}, end coord: {}, meth: {}, aligned seq: {}, refseq: {}'.format(reverse_read.chr,reverse_read.start_coord,reverse_read.end_coord,reverse_read.methylation_string,reverse_read.aligned_ref,reverse_read.seq))

                forward_read_methylation=get_methylation(forward_read,target_objects)
                reverse_read_methylation=get_methylation(reverse_read,target_objects)


                if len(forward_read_methylation)!=0 or len(reverse_read_methylation)!=0: # reads cover at least one target
                    # print(forward_read_methylation)
                    # print(reverse_read_methylation)
                    # print('\n')
                    combined_methylation=forward_read_methylation
                    for cg in reverse_read_methylation:
                        if cg in combined_methylation:
                            if combined_methylation[cg]=='1' and reverse_read_methylation[cg]=='1':
                                combined_methylation[cg]=='1'
                            else:
                                combined_methylation[cg]=='0'
                        else:
                            combined_methylation[cg]=reverse_read_methylation[cg]
                    
                    # populate output line for read
                    new_line=['NA'] * (len(target_objects)+1)
                    new_line[0]=forward_read.read_name

                    for cg in combined_methylation:
                        new_line[target_objects[cg].order]=combined_methylation[cg]
                    
                    output_file.write('\t'.join(new_line)+'\n')

                else: # reads don't cover any targets
                    continue
                    # print('no targets! {}:{}-{}; {}:{}-{}'.format(forward_read.chr,forward_read.start_coord,forward_read.end_coord,reverse_read.chr,reverse_read.start_coord,reverse_read.end_coord))
                
                






    #     # break


    #     # break

    output_file.close()
