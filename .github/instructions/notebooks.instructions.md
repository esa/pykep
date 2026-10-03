---
description: "Use when editing Jupyter notebooks, especially Markdown mathematics rendered by Sphinx, MyST-NB or nbsphinx."
applyTo: "**/*.ipynb"
---

# Notebook conventions

- Follow the repository's Python style and environment instructions for Python code cells and execution.

## Mathematics for Sphinx rendering

The current documentation uses MyST-NB with the `dollarmath` and `amsmath` extensions in `doc/conf.py`, not nbsphinx. Use the following explicit Markdown math boundaries for portable notebook rendering; Jupyter's preview alone is not sufficient validation.

- Use `$...$` for inline mathematics, keeping the expression on one line and without spaces immediately inside the delimiters.
- Use `$$` on separate lines for display mathematics.
- Always leave a blank line before the opening `$$` and after the closing `$$`. Do not attach display equations directly to surrounding prose: a parser can interpret their delimiters as inline math, producing italicized prose, lost spaces and stray dollar signs.
- Keep prose outside math blocks. Use `\text{...}` for words that genuinely belong inside an equation.
- Keep `cases`, `aligned` and matrix environments inside the display delimiters, with matching `\begin{...}` and `\end{...}` and LaTeX `\\` row separators.
- Use notebook Markdown math, not the reStructuredText `.. math::` or `:math:` syntax used in Python docstrings.
- In cell source, use ordinary LaTeX backslashes. When generating notebook JSON, escape backslashes correctly: `\mathbf` in Markdown becomes `\\mathbf` in a JSON string, and a LaTeX `\\` row separator becomes `\\\\`.

### Markdown cell example

```markdown
With `cut` held fixed, the defect Jacobian is

$$
J = \frac{\partial\mathbf{D}}{\partial(\mathbf{X},\mathbf{U},\mathbf{T})}
  = \begin{bmatrix}J_X & J_U & J_T\end{bmatrix}.
$$

Each row corresponds to a defect component; each column corresponds to an input component.
```
