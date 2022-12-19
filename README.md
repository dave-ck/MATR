# MATR

MATR script commissioned by Richard Harris.

### Prerequisites

- Python 3
- Python libraries numpy, pandas, and scipy.

### Behavior

Public data often includes suppressed entries. For the purpose of visualization and analysis, it can be useful to
approximate plausible values for such entries.

Example input table (suppressed values are marked "nan", and would be blank cells in Excel):

```
[[ 138.  495.   62.  695.]
 [  nan  582.   nan 1969.]
 [  nan  659.   nan 1606.]
 [ 943. 1736. 1591. 4270.]]
```

Example guess table -only the cells corresponding to suppressed values will be used, i.e. 12,13,14,15:

```
[[ 138.  495.   62.  695.]
 [  12.  582.   14. 1969.]
 [  13.  659.   15. 1606.]
 [ 943. 1736. 1591. 4270.]]
```

The guess need not be anywhere near the ground truth; in particular, it need not fit constraints from the row and total
columns.

Ground truth table for this example (row and column totals exactly fit):

```
[[ 138.  495.   62.  695.]
 [ 719.  582.  668. 1969.]
 [  86.  659.  861. 1606.]
 [ 943. 1736. 1591. 4270.]]
```

Given an input table and guess table, the script generates an output table such that:

1. All column and row totals are correct for the output.
2. If several solutions would fit (1), then the one which is closest (smallest rectilinear distance) to the guess is
   chosen.

### Directions for use.

1. Ensure prerequisites (above) are correctly installed.
2. Place the script ```matr.py``` in the same directory as ```input.xlsx``` and ```guess.xlsx```.
3. Run the script either:

   a. from the terminal/command prompt, with ```python3 matr.py``` or ```python matr.py``` depending on your
   configuration, or

   b. by double-clicking the script ```matr.py``` in your file explorer.
4. a new file ```output.xlsx``` will appear in the same directory as the script. 
 
