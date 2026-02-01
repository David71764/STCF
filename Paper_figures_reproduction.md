# 论文图像复现说明（用于向作者求助）

本文档对应论文：**Sensitivity study of the charged lepton flavor violating process $\tau \to \gamma\mu$ at STCF**（arXiv:2305.00483）。

我们使用 STCF 全模拟链（KKMC + EvtGen + 探测器模拟 + 重建 + 分析算法）复现论文中的分析流程，并由 `plotting/plot_paper_figs.py` 生成以下图像。当前复现结果与论文仍存在差异，在此按**论文中的图号顺序**列出生成图像及简要说明，便于在飞书或 GitHub 上向作者求助。

**图片所在目录**：`plotting/output/`（运行 `plot_paper_figs.py` 后的默认输出目录）。  
若将本仓库放在 GitHub 上，请保证该目录下存在对应的 PNG 文件，本 Markdown 使用相对路径 `plotting/output/xxx.png`，从仓库根目录打开即可正确显示。

---

## Fig. 1 — 信号侧 $E(\gamma\mu)$ 与 $M(\gamma\mu)$ 二维分布

论文中：框图为信号事例 $E(\gamma\mu)$ 和 $M(\gamma\mu)$ 的二维分布，较大框表示较高密度，红线为信号区（非对称斜椭圆）。

**复现图像（信号样本）：**

![Fig.1 信号样本 E(γμ) vs M(γμ)](plotting/output/EvsM_mugamma_signal.png)

**说明 / 待确认：**  
- 本图由脚本 2D 配置 `EvsM_mugamma`、对信号样本施加 final selection 后绘制，画法为方块大小表示 bin 含量（BOX）。  
- 椭圆信号区由 `ellipse_params.json` 控制，当前若未拟合则图中可能无椭圆线。  
- 若与论文 Fig. 1 在形状、范围或椭圆位置上有差异，请指出可能原因（如选择条件、bin 范围、椭圆参数等）。

---

## Fig. 2 — 运动学分布：(a) $p_\mu$，(b) $E_\gamma$，(c) $\cos\theta_{\gamma\mu}$

论文中：(a) 信号 μ 子动量，(b) 信号光子能量，(c) 信号光子与 μ 子夹角余弦；箭头标记事例选择标准；μ 子动量低端由 MUD 接受度引起。

**复现图像：**

| (a) $p_\mu$ (GeV/c) | (b) $E_\gamma$ (GeV) | (c) $\cos\theta_{\gamma\mu}$ |
|---------------------|----------------------|------------------------------|
| ![pMu](plotting/output/pMu.png) | ![EsigGamma](plotting/output/EsigGamma.png) | ![cosMuGamma](plotting/output/cosMuGamma.png) |

**说明 / 待确认：**  
- 三张图均在“运动学选择、去掉对应变量 cut”的条件下绘制（即 drop 该变量），与论文 Fig. 2 对应。  
- 纵轴为 Events/(0.05 单位)×10³；若与论文纵轴刻度或形状不一致，可能与归一化、bin 宽或样本量有关。

---

## Fig. 3 — $e^+\bar\nu_e\nu_\tau$ 标记模式：(a) $\cos\theta_{\mathrm{sig}\,\gamma,\,\mathrm{tag}}$，(b) $p_{\mathrm{tag}}$，(c) $E_{\mathrm{miss}}$

论文中：$e^+\bar\nu_e\nu_\tau$ 标记下，信号与双 τ 本底的对比；(a) 信号光子与标记带电径迹夹角余弦，(b) 标记带电径迹动量，(c) 丢失能量。

**复现图像：**

| (a) $\cos\theta_{\gamma,\,\mathrm{tag}}$ | (b) $p_{\mathrm{tag}}$ (GeV/c) | (c) $E_{\mathrm{miss}}$ (GeV) |
|------------------------------------------|---------------------------------|--------------------------------|
| ![e_cosSigGamma_TagCharged](plotting/output/e_cosSigGamma_TagCharged.png) | ![e_pTag](plotting/output/e_pTag.png) | ![e_Emiss](plotting/output/e_Emiss.png) |

**说明 / 待确认：**  
- 仅 e-tag（tagMode=1）事例，并施加除图中变量外与论文 Table 1 一致的 cut。  
- (a) 为对数纵轴。若本底/信号相对形状与论文不符，可能与 tag 判定或 E/p 选择有关。

---

## Fig. 4 — $\pi^+\bar\nu_\tau$ 标记模式：(a) $E_{\mathrm{miss}}$，(b) $|\cos\theta_{\mathrm{miss}}|$，(c) $E_{\mathrm{sig}\,\gamma}$，(d) $M^2_{\mathrm{miss}}$

论文中：$\pi^+\bar\nu_\tau$ 标记下，(a)(b)(c) 与双 μ 本底相关，(d) 与双 τ 本底相关。

**复现图像：**

| (a) $E_{\mathrm{miss}}$ | (b) $\lvert\cos\theta_{\mathrm{miss}}\rvert$ | (c) $E_{\mathrm{sig}\,\gamma}$ | (d) $M^2_{\mathrm{miss}}$ |
|-------------------------|-------------------------------------|----------------------------------|----------------------------|
| ![pi_Emiss](plotting/output/pi_Emiss.png) | ![pi_absCosMiss](plotting/output/pi_absCosMiss.png) | ![pi_EsigGamma](plotting/output/pi_EsigGamma.png) | ![pi_M2miss](plotting/output/pi_M2miss.png) |

**说明 / 待确认：**  
- 仅 π-tag（tagMode=2）事例，cut 与 Table 1 一致（除图中被 drop 的变量）。  
- (b) 为对数纵轴。若 $M^2_{\mathrm{miss}}$ 峰位或本底形状与论文有出入，可能与 $\pi^0$ 重建或 missing 四动量定义有关。

---

## Fig. 5 — $\pi^+\pi^0\bar\nu_\tau$ 标记模式：(a) $M^2_{\mathrm{miss}}$，(b) $\cos\theta_{\mathrm{hel}}$，(c) $|\cos\theta_{\mathrm{miss}}|$

论文中：(a)(b) 双 τ 本底，(c) 强子本底中丢失动量方向。

**复现图像：**

| (a) $M^2_{\mathrm{miss}}$ | (b) $\cos\theta_{\mathrm{hel}}$ | (c) $\lvert\cos\theta_{\mathrm{miss}}\rvert$ |
|---------------------------|----------------------------------|-------------------------------------|
| ![rho_M2miss](plotting/output/rho_M2miss.png) | ![rho_cosHel](plotting/output/rho_cosHel.png) | ![rho_absCosMiss](plotting/output/rho_absCosMiss.png) |

**说明 / 待确认：**  
- 仅 ρ-tag（tagMode=3，即 $\pi^+\pi^0\bar\nu_\tau$）事例。  
- 螺旋角在分析中取“信号 τ 静止系中 μ 方向与质心系中 τ 方向”的夹角。若与论文定义或分布形状不一致，请指出。

---

## 2D 图 — 各样本 $E(\gamma\mu)$ vs $M(\gamma\mu)$（供对比）

以下为对**同一 2D 变量**、在施加 final selection 后，按样本分别输出的框图，便于与论文中信号/本底描述对照（论文 Fig. 1 仅给出信号）。

| 信号 Signal | 双 μ Di-μ | 双 τ ττ | 强子 qq̄ |
|-------------|------------|---------|----------|
| ![EvsM_signal](plotting/output/EvsM_mugamma_signal.png) | ![EvsM_dimu](plotting/output/EvsM_mugamma_dimu.png) | ![EvsM_ditau](plotting/output/EvsM_mugamma_ditau.png) | ![EvsM_had](plotting/output/EvsM_mugamma_had.png) |

---

## 目录与文件对应关系（便于 GitHub 检查）

脚本默认将 PNG 输出到 **`plotting/output/`**，文件名与本文档引用一致：

```
plotting/output/
├── EvsM_mugamma_signal.png    # Fig. 1
├── EvsM_mugamma_dimu.png
├── EvsM_mugamma_ditau.png
├── EvsM_mugamma_had.png
├── pMu.png                    # Fig. 2(a)
├── EsigGamma.png              # Fig. 2(b)
├── cosMuGamma.png             # Fig. 2(c)
├── e_cosSigGamma_TagCharged.png   # Fig. 3(a)
├── e_pTag.png                    # Fig. 3(b)
├── e_Emiss.png                   # Fig. 3(c)
├── pi_Emiss.png                  # Fig. 4(a)
├── pi_absCosMiss.png             # Fig. 4(b)
├── pi_EsigGamma.png              # Fig. 4(c)
├── pi_M2miss.png                 # Fig. 4(d)
├── rho_M2miss.png                # Fig. 5(a)
├── rho_cosHel.png                # Fig. 5(b)
└── rho_absCosMiss.png            # Fig. 5(c)
```

生成方式：在项目根目录执行（需先配置 `plotting/samples.json` 中的 root 文件路径）：

```bash
python plotting/plot_paper_figs.py --config plotting/samples.json --outdir plotting/output
```

若将图片放到其他目录（例如 `docs/figures/`），只需修改本 Markdown 中所有 `plotting/output/` 为对应相对路径即可。

---

**致谢**：复现工作基于论文 arXiv:2305.00483 的分析描述与 Table 1 选择条件；若您是该文作者并愿意对上述差异给予指点，将非常感谢。
