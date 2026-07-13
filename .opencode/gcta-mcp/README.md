# GCTA MCP Server

Expose all GCTA functionality to AI tools (Claude, Codex, opencode, etc.) via the MCP (Model Context Protocol) protocol, enabling users to invoke various GCTA analysis functions through natural language conversations.

## Directory Structure

```
GCTA/
├── gcta                  # GCTA binary (Linux x86_64, statically linked)
├── gcta-mcp/
│   ├── server.py         # MCP server entry point
│   ├── requirements.txt  # Python dependencies
│   └── README.md         # This document
├── test.bed / test.bim / test.fam / test.phen  # Sample data
├── MIT_License.txt
└── README.txt
```

`server.py` automatically locates the `../gcta` binary and working directory based on its own location, with no hardcoded paths needed.

## Installation

### 1. Install Python dependencies

```bash
pip3 install mcp
```

### 2. Make the binary executable (if needed)

```bash
chmod +x gcta
```

### 3. Verify

```bash
cd gcta-mcp
python3 server.py
```

If you see the "GCTA MCP Server starting..." message, the server has started successfully.

## Configuration in AI Tools

Replace `<GCTA_DIR>` below with the full path to the `GCTA` folder.

### Claude Desktop

Edit `~/Library/Application Support/Claude/claude_desktop_config.json` (macOS) or the corresponding config file for your platform:

```json
{
  "mcpServers": {
    "gcta": {
      "command": "python3",
      "args": ["<GCTA_DIR>/gcta-mcp/server.py"]
    }
  }
}
```

### opencode

Edit `opencode.json`:

```json
{
  "mcp": {
    "gcta": {
      "type": "local",
      "command": ["python3", "<GCTA_DIR>/gcta-mcp/server.py"],
      "enabled": true
    }
  }
}
```

### Codex (OpenAI)

Edit `~/.codex/config.json`:

```json
{
  "mcpServers": {
    "gcta": {
      "command": "python3",
      "args": ["<GCTA_DIR>/gcta-mcp/server.py"]
    }
  }
}
```

## Environment Variables (Optional)

| Variable | Description | Default |
|----------|-------------|---------|
| `GCTA_BINARY_PATH` | Path to the GCTA binary | `../gcta` (relative to server.py) |
| `GCTA_WORK_DIR` | Working directory (where input/output files reside) | Parent directory of server.py |
| `GCTA_TIMEOUT` | Execution timeout (seconds) | `3600` |

## Tool List

This MCP service provides 23 tools covering all major GCTA analysis functions:

| Tool | Function |
|------|----------|
| `gcta_run` | Generic execution (arbitrary CLI arguments) |
| `gcta_help` | Get GCTA documentation |
| `gcta_check` | Check GCTA configuration |
| `gcta_make_grm` | Build GRM (Genetic Relationship Matrix) |
| `gcta_reml` | REML analysis (estimate variance components/heritability) |
| `gcta_bivariate_reml` | Bivariate REML (estimate genetic correlation) |
| `gcta_hereg` | HE regression analysis |
| `gcta_pca` | Principal component analysis |
| `gcta_mlma` | Mixed linear model association analysis |
| `gcta_cojo` | Conditional and joint analysis |
| `gcta_gsmr` | GSMR Mendelian randomization |
| `gcta_mtcojo` | Multi-trait COJO analysis |
| `gcta_fastgwa` | fastGWA genome-wide association analysis |
| `gcta_fastbat` | fastBAT gene-level association test |
| `gcta_simu_qt` | Simulate quantitative traits |
| `gcta_simu_cc` | Simulate case-control phenotypes |
| `gcta_fst` | Fst population differentiation analysis |
| `gcta_make_bed` | Data management (format conversion, filtering) |
| `gcta_ld_pruning` | LD pruning |
| `gcta_ld_score` | LD score calculation |
| `gcta_acat` | ACAT gene-level test |
| `gcta_list_files` | List files in working directory |
| `gcta_read_file` | Read file contents |

## Usage Examples

Once configured, you can use natural language requests in AI tools:

> "Build a GRM using the test dataset with output prefix test_grm"

The AI will call `gcta_make_grm(bfile="test", out="test_grm")`.

> "Perform REML analysis on test_grm with phenotype file test.phen"

The AI will call `gcta_reml(grm="test_grm", pheno="test.phen", out="test_reml")`.

> "Read the contents of test_reml.hsq"

The AI will call `gcta_read_file(filepath="test_reml.hsq")`.

## Technical Details

- The GCTA binary is in Linux x86_64 ELF format, statically linked (static-pie linked), and runs directly on any Linux x86_64 system
- The MCP service is implemented using the Python `mcp` (FastMCP) library
- All file paths are relative to the working directory (`GCTA_WORK_DIR`, defaulting to the parent directory of server.py)

## License

The GCTA software itself is under the MIT License. This MCP service is free to use and distribute.
