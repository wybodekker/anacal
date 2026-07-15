# anacal

find chemical composition from elemental analysis results

## Properties

|key|value|
|-:|:-|
|  script:|anacal|
|   short:|find chemical composition from elemental analysis results|
|    type:|C|
|  author:|Wybo Dekker|
|   email:|[wybodekker@me.com](mailto:wybodekker@me.com)|
| version:|1.07|
| license:|GNU General Public License|
|   intro:|This script, with the **--install** option, compiles|
|         |and installs a binary anacal in **/usr/local/bin**, with|
|         |documentation in **/usr/local/{man,pdf,html}**.  To this|
|         |end, it creates **anacal.c**, **Makefile**, **README.md** and|
|         |a directory **tests** with input files for testing. These|
|         |can be removed by calling the script with the **--clean**|
|         |option. Any changes in the files will be stored back in|
|         |the script. So if you have applied any changes, run me|
|         |with the **-c** option, and then with the **-I** option to|
|         |install the changes.|

## Options

|option|description|
|:-|:-|
|-I,--install	|install|
|-c,--clean	|clean up after swallowing edited internals and exit|
|-h,--help	|show this help and exit|
|-H,--Help	|show full help and exit|
|-V,--version	|show version and exit|
