---
layout: page
permalink: /cli/
title: Commande line tools
---

# Commande line usage with standalone executables

The command line pipeline for the `DIYABC-RF` inference framework is based on two command-line softwares 
[diyabc](https://github.com/diyabc/diyabc) and [abcranger](https://github.com/diyabc/abcranger).

Please visit the following pages to get the latest release:

- [diyabc](https://github.com/diyabc/diyabc/releases/latest)
- [abcranger](https://github.com/diyabc/abcranger/releases/latest)

## Examples

More details are coming soon.


```sh
diyabc -p ./ -n "'t:4'"
diyabc -p ./ -R "" -m -g 50 -r 100 -t 4
abcranger -t 1000 -j 8 --parameter N1 --chosenscen 1 --noob 50
```
