# 02 — كيمياء النبتة: جرد بدرجات دليل

> المصدر الأساسي: [`../data/compounds.tsv`](../data/compounds.tsv) (**61 مدخلًا**، مُدخلة يدويًّا
> من الأوراق). البنيات الكيميائية في [`../data/01_resolved.csv`](../data/01_resolved.csv)
> (من PubChem) و[`../data/01b_triterpenoids.csv`](../data/01b_triterpenoids.csv) (مبنيّة
> حسابيًّا). الواصفات في [`../data/02_descriptors.csv`](../data/02_descriptors.csv).

---

## 1. الخلاصة في ثلاث جمل

1. **الكيمياء الحقيقية لهذه النبتة هي صابونينات تريتربينية مؤكسدة عند C-23 على هياكل
   30-nor-oleanane وlupane** — لا مشتقّات حمض الأوليانوليك العادية فقط. هذا مثبت بالعزل
   والـ NMR في `salaheldine2019` و`gamal2022`.
2. **كل ما عدا الصابونينات والستيرولات ضعيف الدليل**: 20 مدخلًا من 61 مصدرها مطابقة HPLC/LC-MS
   فقط (`tentative_*`)، و8 مدخلات ملوّثات مخبرية مُعلنة، و5 قلويدات **لم تُعزل من هذا النوع
   أصلًا** بل من أنواع أخرى من الجنس.
3. لذلك: أي حديث عن «المركّبات الفعّالة في الباقل» يجب أن يبدأ من 20 مركّبًا معزولًا، لا من
   قوائم الـ GC-MS.

## 2. توزيع الدليل على الجرد (61 مدخلًا)

| درجة الدليل | العدد | المعنى |
|---|---|---|
| `isolated` | 20 | مُعزل ومُوصَّف بالـ NMR/MS في ورقة أصلية |
| `tentative_hplc` | 20 | مطابقة زمن احتجاز/طيف فقط (أغلبها `aljoufi2022`) |
| `artifact_suspect` | 8 | ملوّثات أو مضافات: BHT، BHA، فثالات، نيتروزومورفولين، anthrone، isopropyl palmitate، **وكافيين** |
| `genus_only` | 5 | قلويدات الجنس (anabasine, lupinine, sparteine, aphylline, aphyllidine) — **لا دليل أوّلي في *A. articulata*** |
| `tentative_lcms` | 3 | مطابقة LC-MS فقط |
| `inferred_aglycone` | 2 | استُنتجت من غلوكوزيدها (β-sitosterol, stigmasterol) |
| `reported` | 2 | مذكورة بلا تفصيل (betaine من LOTUS، 4-acetoxyphenol) |
| `genus_family` | 1 | أوكسالات: موثّقة نسيجيًّا (`doaigey1991`)، بلا كمّية |

## 3. المركّبات المعزولة (العشرون) — جوهر المستودع

### 3.1 الصابونينات والتريتربينات: 30-nor-oleanane وlupane مؤكسدة عند C-23

| المرجع | المركّب (كما كُتب في الورقة) | الهيكل | معرّفنا |
|---|---|---|---|
| `salaheldine2019` | 3β-hydroxy,23-aldehyde-30-norolean-12,20(29)-dien-28-oic acid-28-O-β-D-glucopyranosyl ester (**1**, جديد) | 30-nor-oleanane | `SAP-3` |
| `salaheldine2019` | 3β-O-D-galactopyranosyl-23-aldehyde-30-norolean-12,20(29)-dien-28-oic acid-28-O-β-D-glucopyranosyl ester (**2**, جديد) | 30-nor-oleanane | `SAP-4` |
| `salaheldine2019` | 3β-O-D-xylopyranosyl-30-norolean-12,20(29)-dien-28-oic acid 28-O-β-D-glucopyranosyl ester (**3**, جديد) | 30-nor-oleanane | `SAP-5` |
| `salaheldine2019` | boussingoside E (معروف) | 30-nor-triterpenoid | — (بنيته UNVERIFIED) |
| `gamal2022` | 3β-hydroxy-23-aldehyde-lup-20(29)-ene-28-oic acid (**1**) | lupane | `AGL-5` |
| `gamal2022` | 3β-hydroxy-23-aldehyde-30-nor-olean-12,20(29)-diene-28-oic acid (**2**) | 30-nor-oleanane | `AGL-2` |
| `gamal2022` | 3β-hydroxy-lup-20(29)-ene-23,28-dioic acid (**3**) | lupane | `AGL-6` |
| `gamal2022` | 3β,20α-dihydroxy-30-nor-olean-12-ene-23,28-dioic acid (**4**, جديد) | 30-nor-oleanane | `AGL-4` |
| `gamal2022` | المركّب 3 + 23-O-β-D-glucopyranosyl ester (**5**) | lupane | `SAP-6` |
| `gamal2022` | 3-O-β-D-glucuronopyranosyl-lup-20(29)-ene-23,28-dioic acid 28-O-β-D-glucopyranosyl ester (**6**) | lupane | `SAP-7` |
| `gamal2022` | 3-O-β-D-glucuronopyranosyl-lup-20(29)-ene-23-aldehyde-28-oic acid 28-O-β-D-glucopyranosyl ester (**7**) | lupane | `SAP-8` |
| `metwally2012` | oleanolic acid 3-O-β-D-glucopyranoside | oleanane | `SAP-1` |
| `metwally2012` | 3-O-β-D-glucopyranosyl-28-O-β-D-xylopyranosyl oleanolic acid | oleanane | `SAP-2` |

**النمط الكيميائي الواضح:** C-3 يحمل سكّرًا (Glc / Gal / Xyl / GlcA)، وC-28 إستر غلوكوزي
(بايديسموزيد)، و**C-23 مؤكسد إلى ألدهيد أو حمض** — وهذه آخرها هي السمة غير المعتادة.
هذا النمط يطابق كيمياء العائلة (`mroczek2015`): Amaranthaceae معروفة بـ 30-nor-oleananes
وبتأكسد C-23 (مثل hederagenin وphytolaccagenic acid في الكينوا).

### 3.2 الستيرولات وغيرها

غلوكوزيدات ثلاثة ستيرولات (daucosterol، stigmasterol 3-O-glucoside، sitostanol 3-O-glucoside)
و**proceric acid** (`metwally2012`). الأخير: **لا سجلّ له في PubChem، ولم نعثر على بنيته
المنشورة — مُعلَّم UNVERIFIED في الجرد.** وflavonoids معزولة حديثًا: hyperoside وquercitrin
(`jan2025`).

## 4. البنيات: كيف بُنيت وكيف تحقّقنا منها

PubChem لا يعرف أغلب هذه الصابونينات. لذلك يبنيها
[`../scripts/01b_build_triterpenoids.py`](../scripts/01b_build_triterpenoids.py) من أسلاف
PubChem (حمض الأوليانوليك CID 10494، حمض البيتولينيك CID 64971) **بخطوات معلنة ومُراجَعة**
على ذرّات موسومة (C-3، C-4، C-20، C-23، C-28)، ثم:

1. **تحقّق من الصيغة الجزيئية**: كل بنية من الـ 15 تُقارن بصيغة محسوبة يدويًّا مستقلّة.
   النتيجة: **15/15 مطابقة** (عمود `formula_check`).
2. **تحقّق بصري**: رسم شبكي لكل البنيات في
   [`../data/qc/triterpenoids.png`](../data/qc/triterpenoids.png) — فُحص يدويًّا للتأكّد من
   موضع الألدهيد عند C-23، والمضاعف 20(29)، وموضع السكّريات.
3. **تحقّق خارجي**: `SAP-1` (غلوكوزيد حمض الأوليانوليك) بُني حسابيًّا، وله سجلّ في PubChem
   (CID 165117). التركيب (أول 14 حرفًا من InChIKey) **متطابق تمامًا**:
   `ZNFRITHWVZXJRK`. طبقة الفراغ مختلفة لأن سجلّ PubChem يحمل مراكز فراغية غير محدّدة،
   وبنيتنا محدّدة بالكامل. هذا تحقّق مستقلّ من صحّة طريقة البناء.

### حدود مُعلنة في البناء (عمود `stereo_note`)

- **C-23 مقابل C-24**: C-4 يحمل مثيلين؛ نؤكسد أحدهما. أيّهما «C-23» بالضبط في الـ SMILES
  غير محسوم عندنا. لا أثر على الصيغة ولا الكتلة الدقيقة؛ أثر هامشي على الهندسة ثلاثية الأبعاد.
- **C-20 في `AGL-4`**: التهيئة الفراغية (المنشورة 20α) غير محدّدة في بنيتنا.
- **السكّريات**: الفراغية مأخوذة من أسلاف PubChem (β-D-Glc/Gal/Xyl/GlcA)، وتوصيفات الأوراق
  مقبولة كما كُتبت ولم نُعِد استنتاجها من أطياف NMR (لا نملك الأطياف).

## 5. ماذا تقول الواصفات الدوائية (سكريبت 02)

63 بنية محسوبة. النتيجة محرجة ومفيدة:

- **26 بنية تمرّ فلتري Lipinski وVeber** — لكن انظر من هي: BHT، BHA، فثالات،
  نيتروزومورفولين، وقلويدات *A. aphylla*. أي أن **«المركّبات الشبيهة بالدواء» في هذا الجرد
  أغلبها ملوّثات مخبرية ومركّبات من نوع آخر.** الفلافونويدات (quercetin، kaempferol،
  myricetin) هي الاستثناء الحقيقي.
- **17 صابونينًا/غلايكوزيدًا بكتلة > 550**: كلّها تخرق 2–4 من قواعد Lipinski وحدود Veber.
  أي ادّعاء بأن أحدها «دواء فموي محتمل» يجب أن يشرح أولًا مشكلة الامتصاص.
- **15 بنية ترفع تنبيهات PAINS** (كاتيكولات وبوليفينولات) — وهي بالضبط المركّبات التي تُضيء
  اختبارات DPPH. تفصيل هذا في [`04-evidence-critique.md`](04-evidence-critique.md).

## 6. الفجوة التي تستحقّ العمل

لا يوجد **أي** تحليل HRMS/MS² موجّه للصابونينات في هذه النبتة. كل ما نعرفه عن صابونيناتها
يأتي من أربع عمليات عزل، وكل عملية عزل تلتقط المركّبات الوفيرة فقط. السؤال المفتوح:
**هل الصابونينات المنشورة هي كل ما في النبتة، أم أنها رأس جبل جليد؟**

هذا سؤال يُحسم بحقنة واحدة. القائمة الموجّهة جاهزة:
[`../data/05_ms_inclusion_list.csv`](../data/05_ms_inclusion_list.csv) و
[`../data/05_ms_method_notes.md`](../data/05_ms_method_notes.md).
