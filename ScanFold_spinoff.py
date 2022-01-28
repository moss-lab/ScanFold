### Script to run ScanFold-scan and -fold for webserver
import argparse
import tarfile
import subprocess
import os

def parse_both_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='input filename')
    parser.add_argument('--terminallog', type=str,
                        help='redirect stdout here')
    parser.add_argument('--downloadall', type=str, help='Path and filename write output to')
    parser.add_argument('--name', type=str, default = "UserInput",
                        help='name of data being analyzed')
    parser.add_argument('-t', '--temp', type=str, default="37",
                        help='Folding temperature')
    parser.add_argument('--fasta_index', type=str,
                        help='fasta index file path')
    return parser.parse_known_args()

def parse_scan_args(args):

    #### Parsing arguments ####
    parser = argparse.ArgumentParser()

    parser.add_argument('-s', type=int, default=10,
                        help='step size')
    parser.add_argument('-w', type=int, default=120,
                        help='window size')
    parser.add_argument('-r', type=int, default=30,
                        help='randomizations')
    parser.add_argument('-type', type=str, default='mono',
                        help='randomization type')
    parser.add_argument('-p', type=str, default='off',
                        help='print to screen option (default off:1)')
    parser.add_argument('--print_random', type=str, default='off',
                        help='print to screen option (default off)')
    parser.add_argument('--split', type=str, default = "off",
                        help='name of data being analyzed')

    parser.add_argument('--scan_out_path', type=str,
                        help='ScanFold-Scan output path')
    parser.add_argument('--zscore_wig_file_path', type=str,
                        help='zscore_wig_file_path')
    parser.add_argument('--mfe_wig_file_path', type=str,
                        help='mfe_wig_file_path')
    parser.add_argument('--ed_wig_file_path', type=str,
                        help='ed_wig_file_path')
    parser.add_argument('--pvalue_wig_file_path', type=str,
                        help='pvalue_wig_file_path')
    parser.add_argument('--fasta_file_path', type=str,
                        help='fasta_file path')
    parser.add_argument('--fasta_index', type=str,
                        help='fasta index file path')

    return parser.parse_known_args(args)

def parse_fold_args(args):

    #### Parsing arguments ####
    parser = argparse.ArgumentParser()

    parser.add_argument('-f', type=int, default=-2,
                        help='filter value')
    parser.add_argument('-c', type=int, default=1,
                        help='Competition')
    parser.add_argument('--global_refold', action='store_true',
                        help='global refold oprion')

    ### Required for spinoff ###
    parser.add_argument('--out1', type=str,
                        help='out1 path')
    parser.add_argument('--out2', type=str,
                        help='out3 path')
    parser.add_argument('--out3', type=str,
                        help='out3 path')
    parser.add_argument('--out4', type=str,
                        help='log_file path')
    parser.add_argument('--out5', type=str,
                        help='final_partner_file path')
    parser.add_argument('--out6', type=str,
                        help='bp_track_file path')
    parser.add_argument('--out7', type=str,
                        help='fasta_file path')
    parser.add_argument('--dbn_file_path', type=str,
                        help='dbn_file_path')
    parser.add_argument('--dbn_file_path1', type=str,
                        help='dbn_file_path1')
    parser.add_argument('--dbn_file_path2', type=str,
                        help='dbn_file_path2')
    parser.add_argument('--dbn_file_path3', type=str,
                        help='dbn_file_path3')
    parser.add_argument('--dbn_file_path4', type=str,
                        help='dbn_file_path4')
    parser.add_argument('--structure_extract_file', type=str,
                        help='structure_extract_file path')
    parser.add_argument('--final_partners_wig', type=str,
                        help='final partners wig file path')

    return parser.parse_known_args(args)

if __name__ == "__main__":

    both_args, both_leftover_args = parse_both_args()
    scan_args, scan_leftover_args = parse_scan_args(both_leftover_args)
    fold_args, fold_leftover_args = parse_fold_args(both_leftover_args)

    try:

        scan_args_list = [
            "python", "ScanFold-Scan_spinoff.py",
            '--input', both_args.input,
            '--terminallog', both_args.terminallog,
            '--name', both_args.name,
            '--temp', both_args.temp,
            '--fasta_index', both_args.fasta_index,
        ]
        # the args string to use to call scan is really
        # the args leftover after fold has parsed
        scan_args_list.extend(fold_leftover_args)

        subprocess.run(scan_args_list, check=True, capture_output=True)

        fold_args_list = [
            "python", "ScanFold-Fold_spinoff.py",
            '--input', scan_args.scan_out_path,
            '--terminallog', both_args.terminallog,
            '--name', both_args.name,
            '--temp', both_args.temp,
            '--fasta_index', both_args.fasta_index,
        ]
        # the args string to use to call fold is really
        # the args leftover after scan has parsed
        fold_args_list.extend([scan_leftover_args])

        subprocess.run(fold_args_list, check=True, capture_output=True)

        # generate the tar file of all output

        with open(both_args.downloadall, 'wb') as output_wb:
            tar = tarfile.open(mode="w:gz", fileobj=output_wb)
            tar.add(scan_args.scan_out_path, os.path.basename(scan_args.scan_out_path))
            tar.add(scan_args.zscore_wig_file_path, os.path.basename(scan_args.zscore_wig_file_path))
            tar.add(scan_args.mfe_wig_file_path, os.path.basename(scan_args.mfe_wig_file_path))
            tar.add(scan_args.ed_wig_file_path, os.path.basename(scan_args.ed_wig_file_path))
            tar.add(scan_args.pvalue_wig_file_path, os.path.basename(scan_args.pvalue_wig_file_path))
            tar.add(scan_args.fasta_file_path, os.path.basename(scan_args.fasta_file_path))
            tar.add(scan_args.fasta_index, os.path.basename(scan_args.fasta_index))
            tar.add(fold_args.out1, os.path.basename(scan_args.out1))
            tar.add(fold_args.out2, os.path.basename(scan_args.out2))
            tar.add(fold_args.out3, os.path.basename(scan_args.out3))
            tar.add(fold_args.out4, os.path.basename(scan_args.out4))
            tar.add(fold_args.out5, os.path.basename(scan_args.out5))
            tar.add(fold_args.out6, os.path.basename(scan_args.out6))
            tar.add(fold_args.dbn_file_path4, os.path.basename(scan_args.dbn_file_path4))
            tar.add(fold_args.structure_extract_file, os.path.basename(scan_args.structure_extract_file))
            tar.add(fold_args.final_partners_wig, os.path.basename(scan_args.final_partners_wig))
            tar.close()

    except subprocess.CalledProcessError as e:
        print(e.cmd)
        print(e.stdout)
        print(e.stderr)
        exit(1)