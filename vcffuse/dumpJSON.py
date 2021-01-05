"""
Simple file to pretty-print attributes of a JSON file
"""
import json
import click
import pprint

@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('--jsonfile', '-j', type=str, help='JSON file to parse', required=True)
@click.option('--jkey', '-k', type=str, help='JSON attributes (dict keys) to print out - comma delimited',
              required=False)

def ppjson(jsonfile,jkey):
    # https://stackoverflow.com/a/48111364x
    with open(jsonfile, 'r') as myJSONfile:
        data = myJSONfile.read()
        json_data = json.loads(data)
    if not jkey:
        # we do not have any keys, pretty print the whole file
        pprint.pprint(json_data)
    else:
        # now we have to parse keys:
        for attrkey in jkey.split(","):
            pprint.pprint(json_data[attrkey])

if __name__ == "__main__":
    ppjson()

