#!/usr/bin/env python3

import os
import sys
import pandas as pd


def print_usage():
    print("Usage: %s cicero_annotation_file.txt" % sys.argv[0])


def preprocess_dataframe(df):
    # Order columns (instructions from Scott)
    cols = ["chrA", "posA", "ortA", 
            "chrB", "posB", "ortB",
            "type",
            "geneA", "geneB",
            "featureA", "featureB",
            "matchA", "matchB",
            "repeatA", "repeatB",
            "coverageA", "coverageB",
            "ratioA", "ratioB",
            "readsA", "readsB",
            "qposA", "qposB",
            "total_readsA", "total_readsB",
            "contig"]
    return df[cols]


if __name__ == "__main__":
  
    if len(sys.argv) != 2:
        print_usage()
        sys.exit(1)

    tdt_file = sys.argv[1]
    data = pd.read_table(tdt_file)
    data = preprocess_dataframe(data)

    html_data = data.to_html()
    with open("results.html", "w") as f:
        f.write("<html>\n")
        f.write("<head>\n")
        #f.write('<link rel="stylesheet" type="text/css" href="http://code.jquery.com/ui/1.9.2/themes/base/jquery-ui.css">')
        f.write('<link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.1/css/jquery.dataTables.css">')
        f.write('<link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/responsive/1.0.0/css/dataTables.responsive.css">')
        f.write("</head>\n")
        f.write("<body>\n")
        f.write(html_data+"\n")
        f.write('<script src="https://code.jquery.com/jquery-1.9.1.js"></script>')
        f.write('<script src="http://code.jquery.com/ui/1.9.2/jquery-ui.js"></script>')
        f.write('<script src="https://cdn.datatables.net/1.10.1/js/jquery.dataTables.min.js"></script>')
        f.write('<script src="https://cdn.datatables.net/responsive/1.0.0/js/dataTables.responsive.js"></script>')
        f.write('<script type="text/javascript">\n')
        f.write('$(document).ready(function() {\n')
        f.write('   $(".dataframe").attr("cellspacing", "0").addClass("display").DataTable({"iDisplayLength": 25});\n')
        f.write('});\n')
        f.write("</script>\n")
        f.write("</body>\n")
        f.write("</html>\n")
