# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is Andy Timm's personal website/blog built with [Quarto](https://quarto.org/). The site features technical blog posts primarily focused on statistical methods including variational inference, MRP (Multilevel Regression and Poststratification), BART, and survey methodology.

## Build Commands

**WARNING:** Do not run `quarto render` locally - some posts have Python dependencies that require WSL, and a full render will break the cached outputs. Use `quarto preview` for local development instead.

```bash
# Preview the site locally with live reload (RECOMMENDED)
quarto preview

# Render a single post (if you know it will work locally)
quarto render posts/"Post Name"/filename.qmd
```

Full site rendering happens via GitHub Actions on push to main.

## Site Structure

- `_quarto.yml` - Main site configuration (theme, navbar, output settings)
- `posts/` - Blog posts, each in its own subfolder with a `.qmd` file and related assets
- `posts/_metadata.yml` - Default frontmatter for all posts (freeze: auto, author info, citations)
- `docs/` - Generated output directory (do not edit directly)
- `_freeze/` - Cached computational output from Quarto's freeze feature
- `_extensions/` - Quarto extensions

## Post Structure

Each post lives in `posts/<Post Name>/` containing:
- A `.qmd` file (Quarto markdown) with YAML frontmatter
- Images and other assets referenced in the post
- Posts use `freeze: auto` so computed outputs are cached in `_freeze/`

Post frontmatter typically includes: title, subtitle, date, image, and categories.

## Code Style (R)

When writing R code for posts or analysis:
- Use tidyverse patterns and modern packages
- Use ggplot2 for visualization
- Keep code readable over overly clever
- Document non-obvious choices in markdown or comments
- Consider color blindness in palette choices
- Use assertthat or rlang::abort for input validation

## Technical Notes

- Math rendering uses KaTeX
- The site theme is "journal"
- Posts support both R and Python code execution via Jupyter
- Draft posts can be hidden using `draft: true` in frontmatter (draft-mode: gone)
