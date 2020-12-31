# ADME-T-Prediction
Used for ADME\T attribute prediction
1.ADME\T软件背景介绍
这是一款预测药物分子ADME\T属性的软件，ADME\T属性即药物在生物体内吸收(absorption，A)、分布(distribution，D)、代谢(metabolism，M)、排泄(excretion，E)的规律以及药物的毒性(toxicity，T）。在创新药物研制过程中，药物代谢动力学研究、药效学研究、毒理学研究处于同等重要的地位，已成为药物临床前研究和临床研究的重要组成部分，因此ADME\T属性预测为提高药物研究效率奠定基础。
2.ADME\T各属性及其模型介绍
2.1基本理化性质
LogS：水溶性值的对数。 药物吸收过程的第一步是片剂或胶囊的分解，随后是活性药物的溶解。 显然，低溶解度不利于良好和完全的口服吸收，因此，对该特性的早期测量在药物开发中非常重要。对于LogS属性预测模型我们采用了RF模型，以及2D分子描述符进行构建。
LogD7.4:在pH = 7.4时正辛醇/水分配系数的对数。 为了发挥治疗作用，一种药物必须进入血液循环，然后到达作用部位，合格的药物通常需要在亲脂性和亲水性之间保持平衡，才能溶解在体液中并有效地渗透到生物膜中。因此，估计生理pH下的正辛醇/水分配系数（logD7.4）的值是药物发现早期候选化合物的必备条件。对于LogD7.4属性预测模型我们采用了RF模型，以及2D分子描述符进行构建。
2.2吸收
吸收是指药物从其管理场所进入人体循环系统的过程。对于口服药物来说，肠道是最重要的吸收部位，因此口服药物的人体肠道吸收是必不可少的前提。影响药物吸收的许多因素在不同程度上它们可分为三类：生理因素，如消化系统和循环系统因素；物理化学因素，如离解度和脂溶性；剂型因素，如药物的分解和溶解。
Caco2_Permeability: CaCO-2细胞的通透性。在口服药物到达全身循环之前，它必须通过被动扩散、载体介导的摄取或主动转运过程通过肠道细胞膜。人结肠腺癌细胞株 （Caco-2）作为人肠上皮细胞的一种替代方法，由于其形态和功能相似性，已普遍用于评估体内药物的渗透性。因此，CaCO-2细胞的通透性也是一种合格候选药物的重要指标。对于Caco2_Permeability属性预测模型我们采用了RF模型，以及2D分子描述符进行构建。
Pgp-inhibitor： P-糖蛋白，也称为MDR1或2 ABCB1，是ATP结合盒（ABC）转运蛋白超家族的膜蛋白成员。与hERG通道和CYP3A4一起，它可能是研究最广泛的抗靶标。 实际上，Pgp可能是最混杂的外排转运蛋白，因为它识别许多结构上不同且看似无关的异种生物。 值得注意的是，其中许多也是CYP3A4底物。因此，P-糖蛋白不仅在吸收过程中起着重要作用，而且在分布，代谢和排泄等其他药代动力学过程中也起着重要作用。对于Pgp-inhibitor属性预测模型我们采用了RF模型，以及2D分子描述符进行构建。

HIA： 人体肠道吸收。如前所述，口服药物的人体肠道吸收是其明显疗效的必要前提。而且，口服生物利用度与肠道吸收之间的密切关系也得到了证明，HIA可以在某种程度上被视为口腔生物利用度的替代指标。 为了建立分类模型，定义了正、负化合物。如果化合物A HIA%不到30%，它被标记为阴性(0)；否则，它被标记为阳性(1)。对于HIA属性预测模型我们采用了RF模型，以及MACCS分子描述符进行构建。
F: 人类口服生物利用度。 对于任何通过口服途径给药的药物，口服生物利用度无疑是最重要的药代动力学参数之一，因为它是药物向全身循环传递效率的指标。应用两个阈值（20％和30％）将所有化合物分为阳性和阴性化合物。对于F_20属性预测模型我们采用了RF模型，以及MACCS分子描述符进行构建。对于F_30属性预测模型我们采用了RF模型，以及ECFP6分子描述符进行构建。
2.3分布
通常，药物的分布是血液和组织之间的运输过程。药物从给药部位吸收到血液中后，循环系统将充当转运器，将药物输送到其靶器官，靶组织和靶部位。关于分布的影响因素，主要有药物的理化性质，例如药物的结构特征和脂性，以及人体的生理特征，例如血浆蛋白结合，血流和血管通透性。这些前述因素会导致各种药物的分布差异，并直接影响药物疗效和安全性。
PPB：血浆蛋白结合。众所周知，药物吸收和分配的主要机制之一是通过PPB，因此，药物与血浆中蛋白质的结合对其药代动力学行为有很大影响。一方面，PPB可直接影响口服生物利用度，因为在此过程中，当药物与血清蛋白结合时，药物的自由浓度受到威胁。另一方面，蛋白质-药物复合物可以充当仓库。 因此，有必要在药物开发的早期阶段对其进行评估。对于PPB属性预测模型我们采用了RF模型，以及2D分子描述符进行构建。
VD:发行量。VD是将给药剂量与循环中存在的实际初始浓度联系起来的理论概念，它是描述药物体内分布的重要参数。在实践中，我们可以根据未知化合物的VD值推测其分布特征，例如其结合血浆蛋白的条件，其在体液中的分布量以及在组织中的吸收量。对于VD属性预测模型我们采用了RF模型，以及2D分子描述符进行构建。 
BBB：血脑屏障。BBB是药物的重要药代动力学性质，是其能够或无法穿透血脑屏障。 BBB渗透对于靶向大脑受体的药物很重要。这些药物的例子是抗精神病药，抗癫痫药和抗抑郁药。 对于不针对大脑目标的药物，BBB渗透是不可取的，因为它会导致与CNS相关的不良副作用。对于VD属性预测模型我们采用了SVM模型，以及ECFP2分子描述符进行构建。 
2.4代谢
代谢是生命系统的标志，在其中进行维持体内稳态的复杂生物化学转化。对于所有药物中的约75％，新陈代谢是主要的清除途径之一。通过将新陈代谢系统转化为易于排泄的代谢产物，它已成为抵抗外来有害物质的主要防线。代谢系统高度复杂且适应性强。对于此过程，涉及了许多不同的酶家族，它们通常可以分为两类：微粒体酶，例如对大多数药物重要的细胞色素P450（CYP）酶，以及对少数药物重要的非微粒体酶。因此，在药物开发过程中，对于分子的CYP 450酶底物或抑制剂的识别非常重要。因此，我们研究了几种最常见的与代谢相关的亚型：
CYP1A2抑制剂(CYP_inhibitor_1A2)：采用SVM模型，以及ECFP4分子描述符构建模型
CYP1A2底物(CYP_substrate_1A2):采用RF模型，以及ECFP4分子描述符构建模型
CYP3A4抑制剂(CYP_inhibitor_3A4)：采用SVM模型，以及ECFP4分子描述符构建模型
CYP3A4底物(CYP_substrate_3A4):采用RF模型，以及ECFP4分子描述符构建模型
CYP2C9抑制剂(CYP_inhibitor_2C9)：采用SVM模型，以及ECFP4分子描述符构建模型

CYP2C9底物(CYP_substrate_2C9):采用了RF模型，以及ECFP4分子描述符构建模型
CYP2C19抑制剂(CYP_inhibitor_2C19)：采用SVM模型，以及ECFP2分子描述符构建模型
CYP2C19底物(CYP_substrate_2C19):采用RF模型，以及ECFP4分子描述符构建模型
CYP2D6抑制剂(CYP_inhibitor_2D6)：采用RF模型，以及ECFP2分子描述符构建模型
CYP2D6底物(CYP_substrate_2D6):采了RF模型，以及ECFP2分子描述符构建模型
2.5排泄
对于药物化合物，进入人体后，通常会经历吸收过程，分布过程，代谢过程以及最终的排泄过程。排泄是体内药物或其代谢产物的消除过程。分子的排泄特性会影响药物效率和相应的药物副作用。
CL:药物的清除率。清除率是重要的药代动力学参数，它与分布的体积，半衰期以及给药的频率共同定义。对于CL属性预测模型我们采用了RF模型，以及2D分子描述符进行构建。
T1/2:药物的半衰期。对于T1/2属性预测模型我们采用了RF模型，以及2D分子描述符进行构建。
2.6毒性
hERG: 在心脏去极化和复极化过程中，由hERG编码的电压门控钾通道在调节心脏动作电位和静息电位的交换中起主要作用。hERG阻滞可能导致长时间QT综合征（LQTS），心律不齐和尖端扭转型室速（TdP），导致心悸昏厥甚至猝死。因此，评估与hERG相关的心脏毒性 药物设计/发现流程中的重要一步。对于hERG属性预测模型我们采用了RF模型，以及2D分子描述符进行构建。
HHT:对人体的肝毒性。药物引起的肝损伤是患者安全的高度关注，也是导致药物退出市场的主要原因。 在临床试验中，不良的肝功能常常导致药物开发计划的延迟和昂贵的终止。 因此，尽早发现肝毒性潜力对所有利益相关者都很重要。对于HHT属性预测模型我们采用了RF模型，以及2D分子描述符进行构建。
SkinSen:皮肤敏感性是化学危险性测定和安全性评估的重要毒理学终点。对于SkinSen属性预测模型我们采用了RF模型，以及MACCS分子描述符进行构建。
