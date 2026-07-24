# 接头目录维护与证据工作流

## 设计原则

`config/adapter_catalog.tsv` 保存能够参与 fastp 执行的建库方案；`config/adapter_sequence_reference.tsv` 保存完整PCR引物、flow-cell伪影、反向互补序列及SISPA等协议特异引物。网页、fastp 调度和结果解释读取同一组版本化文件。默认使用 `auto`：fastp 先利用双端重叠剪切，并启用 PE 接头自动识别。只有实验记录明确给出建库试剂盒时，才选择一个手动方案。

不要为了“覆盖更多接头”把所有序列同时传给 `--adapter_fasta`。不相关序列会扩大误匹配和误剪切范围。index、barcode、UMI、完整寡核苷酸结构也不能未经方向换算就当成 read-through trimming sequence。

## 新增或修订一条记录

1. 优先取得建库厂商明确标注的 `Adapter Trimming`、`AdapterRead1` 和 `AdapterRead2` 序列。
2. 若只有软件内置或社区资料，将 `source_level` 标为 `software` 或 `community`，并将 `status` 标为 `review`。
3. 使用稳定的小写 `profile_id`；序列只允许大写 `A/C/G/T`，最短 6 nt。
4. 填写来源 URL、文档版本、核验日期和适用建库名称。不要用搜索结果页作为来源。
5. 执行：

   ```bash
   python3 scripts/helpers/adapter_evidence.py validate \
     --catalog config/adapter_catalog.tsv \
     --reference config/adapter_sequence_reference.tsv
   python3 -m pytest -q
   ```

6. 查看差异并提交：

   ```bash
   git diff -- config/adapter_catalog.tsv
   git add config/adapter_catalog.tsv
   git commit -m "更新接头目录：说明来源和原因"
   ```

## 用户提供接头表的导入审核

用户表应至少提供：建库/试剂盒名称、R1、R2、原始来源或实验记录。逐条完成：

- 去除空格、`5'-`/`-3'` 和修饰符，但保留原始文件不改；
- 检查字符、长度、方向，以及是否误填 index/UMI；
- 与目录做精确和前缀家族比对；
- 厂商来源交叉核验；
- 无法确认方向或用途的项目只进入 `review`，不在网页默认启用。

审核结果和原文件应随同一次 Git 提交保存，从而能够追踪谁在何时、依据什么资料修改了序列。

## 每次运行的可追溯输出

运行目录保存：

```text
.contig_pipeline/runs/<run_id>/
├── parameters.env
├── adapter_catalog.snapshot.tsv
├── adapter_catalog.snapshot.sha256
├── adapter_sequence_reference.snapshot.tsv
└── adapter_sequence_reference.snapshot.sha256
```

每个样本保存：

```text
cleandata/<sample>/fastp_report/
├── <sample>.fastp.html
├── <sample>.fastp.json
├── <sample>.adapter_evidence.tsv
└── <sample>.adapter_reference_scan.tsv
```

证据表区分 `fastp_auto_detection` 与 `configured_fallback_plus_fastp_auto`。后者表示已知接头被配置为双端重叠失败时的后备序列，同时仍启用 PE 自动识别。序列匹配只能说明它与某接头家族一致，不能单凭 FASTQ 反推出唯一试剂盒。

参考扫描默认抽取每个 mate 的前 100,000 条 reads，对版本化参考序列做精确匹配并分别记录 5′、3′和任意位置命中。该扫描用于发现协议特异标签和伪影，不会把 `reference_only`、`protocol_specific` 或 `do_not_use` 条目自动传给 fastp。

## 发布、校正与回退

稳定更新通过测试后创建提交和版本标签，再同步到正式目录。发布后发现序列错误时，不修改旧运行快照；修订目录、增加说明并创建新提交。需要恢复软件版本时用 `git revert <commit>` 创建一条可审计的反向提交，不使用 `git reset --hard` 改写历史。
