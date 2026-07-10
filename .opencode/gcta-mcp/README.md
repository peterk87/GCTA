# GCTA MCP Server

通过 MCP (Model Context Protocol) 协议，将 GCTA 的所有功能暴露给 AI 工具（Claude、Codex、opencode 等），使用户能够通过自然语言对话调用 GCTA 软件的各种分析功能。

## 目录结构

```
GCTA/
├── gcta                  # GCTA 二进制（Linux x86_64, 静态链接）
├── gcta-mcp/
│   ├── server.py         # MCP 服务入口
│   ├── requirements.txt  # Python 依赖
│   └── README.md         # 本文档
├── test.bed / test.bim / test.fam / test.phen  # 示例数据
├── MIT_License.txt
└── README.txt
```

`server.py` 通过自身位置自动定位 `../gcta` 二进制和工作目录，无需任何硬编码路径。

## 安装

### 1. 安装 Python 依赖

```bash
pip3 install mcp
```

### 2. 赋予二进制可执行权限（如果需要）

```bash
chmod +x gcta
```

### 3. 验证

```bash
cd gcta-mcp
python3 server.py
```

如果看到 "GCTA MCP Server starting..." 消息，说明服务器启动成功。

## 在 AI 工具中配置

将下方的 `<GCTA_DIR>` 替换为 `GCTA` 文件夹的完整路径。

### Claude Desktop

编辑 `~/Library/Application Support/Claude/claude_desktop_config.json`（macOS）或对应平台的配置文件：

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

编辑 `opencode.json`：

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

编辑 `~/.codex/config.json`：

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

## 环境变量（可选）

| 变量 | 说明 | 默认值 |
|------|------|--------|
| `GCTA_BINARY_PATH` | GCTA 二进制文件路径 | `../gcta`（相对于 server.py） |
| `GCTA_WORK_DIR` | 工作目录（输入/输出文件所在） | server.py 的上级目录 |
| `GCTA_TIMEOUT` | 执行超时时间（秒） | `3600` |

## 工具列表

本 MCP 服务提供以下 23 个工具，覆盖 GCTA 的全部主要分析功能：

| 工具 | 功能 |
|------|------|
| `gcta_run` | 通用执行（任意命令行参数） |
| `gcta_help` | 获取 GCTA 文档 |
| `gcta_check` | 检查 GCTA 配置 |
| `gcta_make_grm` | 构建 GRM（遗传关系矩阵） |
| `gcta_reml` | REML 分析（估计方差组分/遗传力） |
| `gcta_bivariate_reml` | 双变量 REML（估计遗传相关性） |
| `gcta_hereg` | HE 回归分析 |
| `gcta_pca` | 主成分分析 |
| `gcta_mlma` | 混合线性模型关联分析 |
| `gcta_cojo` | 条件和联合分析 |
| `gcta_gsmr` | GSMR 孟德尔随机化 |
| `gcta_mtcojo` | 多性状 COJO 分析 |
| `gcta_fastgwa` | fastGWA 全基因组关联分析 |
| `gcta_fastbat` | fastBAT 基因水平关联检验 |
| `gcta_simu_qt` | 模拟数量性状 |
| `gcta_simu_cc` | 模拟病例-对照表型 |
| `gcta_fst` | Fst 群体分化分析 |
| `gcta_make_bed` | 数据管理（格式转换、过滤） |
| `gcta_ld_pruning` | LD 剪枝 |
| `gcta_ld_score` | LD 评分计算 |
| `gcta_acat` | ACAT 基因水平检验 |
| `gcta_list_files` | 列出工作目录文件 |
| `gcta_read_file` | 读取文件内容 |

## 使用示例

配置完成后，在 AI 工具中可以直接用自然语言请求：

> "使用 test 数据集构建 GRM，输出前缀为 test_grm"

AI 会调用 `gcta_make_grm(bfile="test", out="test_grm")`。

> "对 test_grm 进行 REML 分析，表型文件为 test.phen"

AI 会调用 `gcta_reml(grm="test_grm", pheno="test.phen", out="test_reml")`。

> "读取 test_reml.hsq 文件的内容"

AI 会调用 `gcta_read_file(filepath="test_reml.hsq")`。

## 技术细节

- GCTA 二进制为 Linux x86_64 ELF 格式，静态链接（static-pie linked），可在任何 Linux x86_64 系统上直接运行
- MCP 服务使用 Python `mcp` (FastMCP) 库实现
- 所有文件路径相对于工作目录（`GCTA_WORK_DIR`，默认为 server.py 的上级目录）

## 许可证

GCTA 软件本身遵循 MIT 许可证。本 MCP 服务可自由使用和分发。
