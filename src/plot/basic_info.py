# Dataset Metrics x Interaction Metrics (Jaccard, Simpson) について以下の三つのグラフをそれぞれ作成する
# 1. 上位 0.25 Quantileに対するクロス集計図のヒートマップ
# 2. hueを上位0.25 Dataset Metricsとした時のInteraction Metricsのバイオリンプロット
# 3. 0.25, 0.5, 0.75 Quantileに対するFisherの正確検定の結果をヒートマップにまとめる

# 入力条件
# HepG2 Over 1000 interaction
# HepG2 Under 1000 interaction
# K562 Over 1000 interaction
# K562 Under 1000 interaction


if __name__ == "__main__":
    from src.eclip.sampleset import SampleSetECLIP
    from src.eclip.encodecondition import ECLIP_SAMPLESETS

    for condition in ECLIP_SAMPLESETS:
        sampleset = SampleSetECLIP(condition)
        print(sampleset)
        print(sampleset.report.shape)
