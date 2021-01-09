class SnpEffParser:
    """
        We are processing an snpEff annotation entry (not a VCF line, but the ANN= stuff)
        and returning with a list of dictionaries
    """

    def __init__(self):
        self._ann_list = list()  # list of annotations
        self._ann_keys = list()

    @property
    def ann_keys(self):
        return self._ann_keys

    @property
    def ann_list(self):
        return self._ann_list

    def add_ann_list(self, ann: dict):
        self._ann_list.append(ann)

    def add_ann_keys(self, ann_line: str):
        """
        Parsing the ##INFO=<ID=ANN line
        """
        entry_list = ann_line.split("'")[1]
        entries = entry_list.split("|")
        for item in entries:
            self._ann_keys.append(item.strip())

    def parse_se_ann(self, line, regexp):
        """
        Annotation entries are delimited by commas, process each entry, and add to a list if contains gene_fusion
        """
        ann_dict = {}
        entries = line.split("ANN=")[1]
        entries = entries.split(",")
        for tr_ann in entries:  # get the transcript items
            if regexp.match(tr_ann):
                # chop up values into a list
                ann_list = tr_ann.split("|")
                # add values into a dict
                for i in range(0, len(self.ann_keys)):
                    ann_dict[self.ann_keys[i]] = ann_list[i]
                # add dict to collection (list)
                self._ann_list.append(ann_dict)
