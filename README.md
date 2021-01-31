# anacal
|     key | description
|     ---:|:---
|  script | anacal - find chemical composition from elemental analysis results
|    type | bash
|  author | Wybo Dekker
|   email | wybo@dekkerdocumenten.nl
| version | 0.00
| license | GNU General Public License

Given, for a chemical composition:
- the elemental percentages for (some of the) elements and
- a set of bruto formulae, one for each possible component
anacal calculates compositions with one or more of the formulas
resulting in a good (that is: better than a certain standard deviation) fit
between calculated and experimental element percentages.
