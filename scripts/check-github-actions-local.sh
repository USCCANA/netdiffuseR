#!/usr/bin/env bash
#
# check-github-actions-local.sh - Simulate GitHub Actions checks locally
#
# This script replicates the GitHub Actions workflow locally using Docker
# to test changes before pushing to GitHub.
#
# Usage:
#   ./scripts/check-github-actions-local.sh [ubuntu|macos|windows|all]
#
# Requirements:
#   - Docker installed and running
#   - Run from repository root

set -e
set -u

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Configuration
PLATFORM="${1:-ubuntu}"
PACKAGE_NAME=$(awk '/^Package:/ {print $2}' DESCRIPTION)
PACKAGE_VERSION=$(awk '/^Version:/ {print $2}' DESCRIPTION)

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}GitHub Actions Local Check${NC}"
echo -e "${GREEN}========================================${NC}"
echo "Package:  ${PACKAGE_NAME} ${PACKAGE_VERSION}"
echo "Platform: ${PLATFORM}"
echo ""

case "$PLATFORM" in
  ubuntu|linux)
    echo -e "${YELLOW}[1/3] Testing Ubuntu (r-release)...${NC}"
    
    # Use R-hub Ubuntu container (simulates ubuntu-latest with R release)
    docker run --rm \
      -v "$(pwd):/workspace" \
      -w /workspace \
      ghcr.io/r-hub/containers/ubuntu-release:latest \
      bash -c "
        set -e
        echo '>>> Installing system dependencies...'
        apt-get update > /dev/null 2>&1
        apt-get install -y libssl-dev libcurl4-openssl-dev libglpk-dev libgmp3-dev libxml2-dev > /dev/null 2>&1
        
        echo '>>> Installing R package dependencies...'
        Rscript -e \"install.packages(c('Rcpp', 'RcppArmadillo', 'sna', 'network', 'networkDynamic', 'Matrix', 'MASS', 'MatchIt', 'SparseM', 'boot', 'igraph', 'viridisLite', 'rcmdcheck', 'remotes', 'knitr'), repos='https://packagemanager.posit.co/cran/__linux__/noble/latest', quiet=TRUE)\"
        
        echo '>>> Running R CMD check...'
        Rscript -e \"rcmdcheck::rcmdcheck(args = c('--no-manual', '--as-cran'), error_on = 'warning', check_dir = 'check')\"
      " && echo -e "${GREEN}✓ Ubuntu check passed${NC}" || echo -e "${RED}✗ Ubuntu check failed${NC}"
    ;;
    
  macos|mac)
    echo -e "${YELLOW}macOS checks require macOS host - skipped${NC}"
    echo "To test on macOS, run: R CMD build . && R CMD check --as-cran *.tar.gz"
    ;;
    
  windows|win)
    echo -e "${YELLOW}Windows checks require Windows host - skipped${NC}"
    echo "To test on Windows, use rhub or GitHub Actions"
    ;;
    
  all)
    echo -e "${BLUE}Running all available checks...${NC}"
    $0 ubuntu
    ;;
    
  *)
    echo -e "${RED}Unknown platform: $PLATFORM${NC}"
    echo "Usage: $0 [ubuntu|macos|windows|all]"
    exit 1
    ;;
esac

echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}Local check completed!${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo "Note: This simulates GitHub Actions but may have minor differences."
echo "For exact GitHub Actions behavior, push to a branch and create a PR."
echo ""
