# Copilot Instructions for netdiffuseR

## Documentation Style

When writing roxygen2 documentation, prefer markdown syntax over LaTeX:

- Use `` `code` `` instead of `\code{code}`
- Use `` `[function()]` `` instead of `\code{\link{function}}`
- Use markdown bullet lists instead of `\itemize{}`
- Use markdown formatting for emphasis and structure

## Code Style

- Follow existing code patterns and conventions in the package
- Use meaningful variable names and keep functions focused
- Add comments only when they match the existing style or explain complex logic
- Prefer using existing libraries over adding new dependencies

## Testing

- Add comprehensive tests for new functions
- Include edge cases and error conditions
- Test with the existing sample datasets when possible
- Validate input parameters and handle errors gracefully

## Examples

- Use `\dontrun{}` only when examples take a long time to run
- Prefer examples that can execute quickly for CRAN checks
- Use existing package datasets in examples when available