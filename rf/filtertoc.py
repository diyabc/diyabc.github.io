"""
Set all headers to level 1
"""

import panflute as pf
import sys
import re


def action(elem, doc):
    if isinstance(elem, pf.Div):
        if 'custom-style' in elem.attributes:
            if 'toc' in elem.attributes['custom-style']:
                del elem.attributes['custom-style']
                return []
    if isinstance(elem, pf.Header):
        sys.stderr.write(elem.content[0].text + "\n")
        if not re.match("^(\d+\.)+", elem.content[0].text):
            elem.content.append(pf.RawInline('\n{: .no_toc }'))
            return elem
        # else:
        #     del elem.content[0]
        #     return elem
            # elem.content.append(pf.Str("\n{:toc}"))
            # sys.stderr.write("replaced" + "\n")
            # return pf.Plain(*elem.content)
        # sys.stderr.write(elem.identifier + "\n")
        # sys.stderr.write(str(elem.to_json()) + "\n\n")

        # elem.identifier = elem.content[0].text.replace(
        #     '.', '') + "-" + elem.identifier


if __name__ == '__main__':
    pf.run_filter(action)
