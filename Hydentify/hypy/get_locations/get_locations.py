'''
get_locations.py: Run SpaceEx using hypy and print out the reachable locations to a file

Stanley Bak
9-2015
'''

HYPY_FOLDER = './hypy/hybridpy'

SPACEEX_SCRIPT = './spaceex/spaceex'

import sys
import os

sys.path.append(HYPY_FOLDER)
import hybridpy as hypy

os.environ["SPACEEX_BIN"] = os.path.realpath(SPACEEX_SCRIPT)

def main():
    '''runs spaceex on model and outputs result'''

    if len(sys.argv) != 3:
        print "Expected 2 params (got " + str(len(sys.argv)-1) + "): [.xml model file] [output path]"
        exit(1)
    
    model = sys.argv[1]
    output_path = sys.argv[2]

    e = hypy.Engine()
    e.set_model(model)
    e.set_tool('spaceex')
    e.set_create_result(True)

    # uncomment this to print the tool's output to stdout
    #e.set_print_terminal_output(True)

    code = e.run(run_hyst=False, make_image=False)

    if code != hypy.RUN_CODES.SUCCESS:
        print 'hypy error: ' + str(code)
        exit(1)

    result = e.get_result()
    locs = result['locations'].keys()
    fixpoint = result.get('fixpoint', False)
    output_separator = '\n'

    with open(output_path, 'w') as f:
        if fixpoint == True:
            f.write('fixpoint found\n')
        else:
            f.write('no fixpoint found\n')
        
        f.write(output_separator.join(locs))


main()
